#!/usr/bin/env python3
"""nf-core/deepmutscan - mutational-signature fit of the sequencing/library error spectrum, one FDM pool.

Spectra (96 pyrimidine-centred trinucleotide channels):
  fdm     false-double error rate of the pooled libraries:  sum F / sum M per channel
  wt      wildtype-template error rate (if the sample has a wildtype library): sum counts / sum coverage
  excess  fdm minus wt, after removing genuine double mutants from F (see below)
Rates are per site and read, i.e. per opportunity. They are put on GRCh37 trinucleotide composition (rate x
GRCh37 count) so that they can be compared with COSMIC-style signatures, which are defined on that composition.

Genuine double mutants. A false-double read is an error only if the second variant is not really on the
molecule. The evidence step measures how often two programmed variants co-occur on reads where both are
multi-nucleotide (so error-free): kappa = genuine-double reads / (library frequency x anchor exposure).
E[genuine-double reads in F_j] = kappa * f_j * M_j, with f_j the SNV's library frequency minus the channel's
FDM error rate. These reads are removed before the excess is formed.

Two-stage fit on each spectrum, against a candidate signature set (default: COSMIC v3.4 library/sequencing
artefacts plus DNA-damage chemistry that can act in vitro; clock-like signatures excluded):
  1. screen: non-negative least squares (NNLS) on all candidates, repeated on bootstrap resamples of the
     sites within each channel; a candidate passes if its exposure is >= min_exposure in >= min_freq of them.
  2. refit: forward selection among the passing candidates by reconstruction cosine, skipping any candidate
     whose cosine to an already chosen signature is >= max_cosine; the number kept is the elbow (kneedle) of
     the cosine curve. Exposures of the final set get bootstrap 95 % intervals.
Also reported: the share of reads the reconstruction places in the right channel (1 - total-variation
distance, in the libraries' own composition).

Nextflow template: values below are interpolated; the Python avoids dollar signs.
"""
import json, platform
import numpy as np
import pandas as pd

POOL = "${meta.id}"
EVIDENCE = "${evidence}".split()
SUMMARIES = "${summaries}".split()
WT_EVIDENCE = "${wt_evidence}".split()
RATES = "${rates}".split()
DESCRIPTIONS = "${descriptions}"
CANDIDATES = "${candidates}"
TRINUC = "${trinucleotides}"
MAX_COSINE = float("${max_cosine}")
BOOTSTRAPS = int("${bootstraps}")
GENUINE = "${genuine_doubles}".lower() == "true"
MIN_EXPOSURE, MIN_FREQ, SEED = 0.02, 0.75, 1


def nnls(A, b, tol=None):
    """Lawson-Hanson non-negative least squares, min ||A x - b|| subject to x >= 0 (numpy only).
    Least-squares subproblems are solved on the columns of A (not the normal equations), for conditioning."""
    m, n = A.shape
    if tol is None:
        tol = 10 * max(m, n) * np.finfo(float).eps * np.linalg.norm(A, 1) * max(1.0, np.abs(b).max())
    x, P = np.zeros(n), np.zeros(n, bool)
    w = A.T @ (b - A @ x)
    for _ in range(30 * n + 30):
        cand = (~P) & (w > tol)
        if not cand.any():
            break
        j = np.argmax(np.where(cand, w, -np.inf))
        P[j] = True
        while True:
            idx = np.flatnonzero(P)
            z = np.zeros(n)
            z[idx] = np.linalg.lstsq(A[:, idx], b, rcond=None)[0]
            if (z[idx] > 0).all():
                x = z
                break
            neg = idx[z[idx] <= 0]
            alpha = np.min(x[neg] / (x[neg] - z[neg]))
            x = x + alpha * (z - x)
            P &= x > tol
            x[~P] = 0
            if not P.any():
                break
        w = A.T @ (b - A @ x)
    return x


def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / na / nb) if na > 0 and nb > 0 else 0.0


## ---- reference data -----------------------------------------------------------------------------
SIG = pd.read_csv(CANDIDATES, index_col=0)
CH = list(SIG.index)                                        # 96 channels, e.g. A[C>A]A
SIG = SIG / SIG.sum()
tri = pd.read_csv(TRINUC, index_col=0).iloc[:, 0]
G37 = np.array([tri[c[0] + c[2] + c[6]] for c in CH], float)
CI = {c: i for i, c in enumerate(CH)}

## ---- site table -----------------------------------------------------------------------------------
ev = pd.concat([pd.read_csv(f) for f in EVIDENCE], ignore_index=True)
ev = ev[ev["context"].notna() & np.isfinite(ev["M"]) & np.isfinite(ev["F"])]
S = ev.groupby("base_mut").agg(context=("context", "first"), codon=("codon", "first"), F=("F", "sum"),
                                M=("M", "sum"), counts=("counts", "sum"), cov=("cov", "sum")).reset_index()
S = S[(S["M"] > 0) & S["context"].isin(CI)].reset_index(drop=True)
S["ch"] = S["context"].map(CI).astype(int)
S["lib_freq"] = S["counts"] / S["cov"]

kappa, gd_reads, gd_exposure = 0.0, 0.0, 0.0
for f in SUMMARIES:
    sm = pd.read_csv(f)
    if "gd_reads" in sm and "gd_exposure" in sm:
        gd_reads += float(np.nansum(sm["gd_reads"]))
        gd_exposure += float(np.nansum(sm["gd_exposure"]))
if GENUINE and gd_exposure > 0:
    kappa = gd_reads / gd_exposure
chsum = lambda v: np.bincount(S["ch"], weights=v, minlength=96)
r_fdm = chsum(S["F"]) / np.maximum(chsum(S["M"]), 1e-300)
S["f_true"] = np.clip(S["lib_freq"] - r_fdm[S["ch"]], 0, None)
# DMSANALYSIS_FDM_ESTIMATE reports G per variant: the genuine-double reads it removed, with every anchor
# weighted by kappa at its own codon distance. Prefer it, so that the spectrum is built on exactly the
# evidence the correction used; the flat kappa above is the fallback when the rates table has no G column
# (--subtract_genuine_doubles false, or a rates table from an older run).
G = None
if GENUINE and RATES:
    rt = pd.read_csv(RATES[0])
    if "G" in rt and float(np.nansum(rt["G"])) > 0:
        G = pd.Series(rt["G"].values, index=rt["base_mut"]).reindex(S["base_mut"]).fillna(0).values
S["F_gd"] = kappa * S["f_true"] * S["M"] if G is None else G
S["F_adj"] = np.clip(S["F"] - S["F_gd"], 0, None)
GD_SOURCE = "flat kappa" if G is None else "kappa(d), from the pooled fit"

def _gd_summary(S, applied, source, kappa, gd_reads):
    """The genuine-double subtraction, library-wide and as a range over the channels. A channel needs at
    least 10 false-double reads to be quoted, so the range is not set by one or two reads."""
    f = np.bincount(S["ch"], weights=S["F"], minlength=96)
    g = np.bincount(S["ch"], weights=S["F_gd"], minlength=96)
    ok = f >= 10
    sh = (g[ok] / f[ok]) if ok.any() else np.array([])
    return {"applied": applied, "source": source, "kappa": round(kappa, 5), "reads": gd_reads,
            "share_of_F": round(float(S["F_gd"].sum() / max(S["F"].sum(), 1e-300)), 4),
            "share_lo": round(float(sh.min()), 4) if len(sh) else None,
            "share_hi": round(float(sh.max()), 4) if len(sh) else None,
            "share_n_channels": int(ok.sum())}


HAS_WT = len(WT_EVIDENCE) > 0
if HAS_WT:
    wt = pd.concat([pd.read_csv(f) for f in WT_EVIDENCE], ignore_index=True)
    w = wt.groupby("base_mut").agg(c=("counts", "sum"), v=("cov", "sum"))
    S["wt_count"] = w["c"].reindex(S["base_mut"]).fillna(0).values
    # variants not observed in the wildtype library still had coverage: take their codon's (or the nearest
    # codon's) median wildtype coverage
    cod = wt.groupby("codon")["cov"].median()
    if len(cod):
        cc = np.array(cod.index, float)
        near = cod.values[np.abs(S["codon"].values[:, None] - cc[None, :]).argmin(axis=1)]
        S["wt_cov"] = w["v"].reindex(S["base_mut"]).values
        S["wt_cov"] = S["wt_cov"].fillna(pd.Series(near, index=S.index))
    else:
        HAS_WT = False
if HAS_WT:
    S = S[S["wt_cov"] > 0].reset_index(drop=True)
    r_wt = np.bincount(S["ch"], weights=S["wt_count"], minlength=96) / np.maximum(np.bincount(S["ch"], weights=S["wt_cov"], minlength=96), 1e-300)
    S["X"] = S["F_adj"] - S["M"] * r_wt[S["ch"]]

## the false-double spectrum is built on F_adj, the same genuine-double-subtracted reads the correction
## uses, so that the panel, the tile above it and the rates applied downstream are one quantity
KINDS = {"fdm": ("F_adj", "M")}
if HAS_WT:
    KINDS.update({"wt": ("wt_count", "wt_cov"), "excess": ("X", "M")})
GROUPS = [np.flatnonzero(S["ch"].values == k) for k in range(96)]
GROUPS = [g for g in GROUPS if len(g)]
rng = np.random.default_rng(SEED)
BOOT = [np.concatenate([rng.choice(g, len(g)) for g in GROUPS]) for _ in range(BOOTSTRAPS)]


def spectrum(num, den, idx=None):
    n, d = (num, den) if idx is None else (num[idx], den[idx])
    ch = S["ch"].values if idx is None else S["ch"].values[idx]
    r = np.bincount(ch, weights=n, minlength=96) / np.maximum(np.bincount(ch, weights=d, minlength=96), 1e-300)
    x = np.clip(r, 0, None) * G37
    return x / x.sum() if x.sum() > 0 else x


def fit(x, sigs):
    A = SIG[sigs].values
    w = nnls(A, x)
    return w, cosine(A @ w, x)


def two_stage(num, den):
    x = spectrum(num, den)
    cands = list(SIG.columns)
    W = []
    for idx in BOOT:
        w, _ = fit(spectrum(num, den, idx), cands)
        W.append(w / w.sum() if w.sum() > 0 else w)
    W = np.array(W) if W else np.zeros((0, len(cands)))
    w0, _ = fit(x, cands)
    med = np.median(W, axis=0) if len(W) else w0 / max(w0.sum(), 1e-300)
    freq = (W >= MIN_EXPOSURE).mean(axis=0) if len(W) else (w0 / max(w0.sum(), 1e-300) >= MIN_EXPOSURE).astype(float)
    screen = [{"signature": s, "median": round(float(m), 4), "freq": round(float(fq), 3),
               "passes": bool(fq >= MIN_FREQ and m >= MIN_EXPOSURE)} for s, m, fq in zip(cands, med, freq)]
    screen.sort(key=lambda r: -r["median"])
    left = [r["signature"] for r in screen if r["passes"]]
    chosen, curve, skipped = [], [], []
    while left:
        best = max(left, key=lambda s: fit(x, chosen + [s])[1])
        left.remove(best)
        rho = max([cosine(SIG[best].values, SIG[c].values) for c in chosen], default=0.0)
        if rho >= MAX_COSINE:
            skipped.append({"signature": best, "max_cosine": round(rho, 3)})
            continue
        chosen.append(best)
        curve.append(round(fit(x, chosen)[1], 4))
    if len(curve) <= 2:
        elbow = len(curve)
    else:
        k = np.arange(1, len(curve) + 1)
        kn = (k - 1) / (k[-1] - 1)
        cn = (np.array(curve) - curve[0]) / max(curve[-1] - curve[0], 1e-12)
        elbow = int(k[np.argmax(cn - kn)])
    final = chosen[:elbow]
    # the burden behind the normalised shape: the measured rate per channel, before the GRCh37
    # rescaling, in reads per million bases of exposure. Reporting only - the report shows it on hover.
    nch = np.bincount(S["ch"], weights=num, minlength=96)
    dch = np.bincount(S["ch"], weights=den, minlength=96)
    res = {"observed": [round(float(v), 6) for v in x],
           "rate_per_mb": [round(float(1e6 * n / d), 3) if d > 0 else None for n, d in zip(nch, dch)],
           "burden_per_mb": round(float(1e6 * nch.sum() / max(dch.sum(), 1e-300)), 2),
           "screen": screen, "order": chosen, "curve": curve,
           "skipped_collinear": skipped, "elbow": elbow, "signatures": final, "exposures": {}, "contributions": {},
           "cosine": None, "reads_accounted": None}
    if not final:
        return res
    w, c = fit(x, final)
    res["cosine"] = round(c, 4)
    tot = w.sum()
    Wb = np.array([fit(spectrum(num, den, idx), final)[0] for idx in BOOT]) if BOOT else np.zeros((0, len(final)))
    Wb = Wb / np.maximum(Wb.sum(axis=1, keepdims=True), 1e-300) if len(Wb) else Wb
    for i, s in enumerate(final):
        res["exposures"][s] = {"w": round(float(w[i] / tot), 4),
                               "lo": round(float(np.quantile(Wb[:, i], .025)), 4) if len(Wb) else None,
                               "hi": round(float(np.quantile(Wb[:, i], .975)), 4) if len(Wb) else None}
        res["contributions"][s] = [round(float(v), 6) for v in SIG[s].values * w[i] / tot]
    # reads accounted for: observed reads per channel vs the reconstruction moved back to the libraries'
    # own composition (divide by GRCh37 opportunity, multiply by the channel's exposure)
    rec = SIG[final].values @ w
    obs = np.clip(np.bincount(S["ch"], weights=num, minlength=96), 0, None)
    opp = np.bincount(S["ch"], weights=den, minlength=96)
    exp_ = rec / G37 * opp
    if obs.sum() > 0 and exp_.sum() > 0:
        res["reads_accounted"] = round(float(np.minimum(obs / obs.sum(), exp_ / exp_.sum()).sum()), 4)
    return res


out = {"pool": POOL, "channels": CH, "n_sites": int(len(S)), "n_libraries": len(EVIDENCE),
       "has_wildtype": HAS_WT, "candidates": list(SIG.columns), "max_cosine": MAX_COSINE, "bootstraps": BOOTSTRAPS,
       "genuine_doubles": _gd_summary(S, GENUINE, GD_SOURCE, kappa, gd_reads),
       "spectra": {}}
if HAS_WT:
    out["excess_share_of_F"] = round(float(S["X"].sum() / max(S["F"].sum(), 1e-300)), 4)
for kind, (a, b) in KINDS.items():
    out["spectra"][kind] = two_stage(S[a].values.astype(float), S[b].values.astype(float))

## one-line glossary for every candidate signature, fitted or not: the report lists the whole set it
## screened, so a reader can see what was on offer as well as what was chosen
notes = {}
try:
    nd = pd.read_csv(DESCRIPTIONS)
    used = sorted(set(SIG.columns)
                  | {sg for r in out["spectra"].values() for sg in r.get("final", [])}
                  | {sg for r in out["spectra"].values() for sg in r.get("order", [])})
    for _, row in nd[nd["signature"].isin(used)].iterrows():
        notes[row["signature"]] = {"label": row["label"], "description": row["description"]}
except Exception as exc:                                        # a user-supplied candidate set need not have notes
    print("no signature descriptions:", exc)
out["descriptions"] = notes

with open("error_signatures.json", "w") as fh:
    json.dump(out, fh, indent=1)
rows = {"channel": CH}
for kind, r in out["spectra"].items():
    rows[f"{kind}_observed"] = r["observed"]
    for s, v in r["contributions"].items():
        rows[f"{kind}_{s}"] = v
pd.DataFrame(rows).to_csv("error_spectra.csv", index=False)
fit_rows = []
for kind, r in out["spectra"].items():
    for s, e in r["exposures"].items():
        fit_rows.append({"pool": POOL, "spectrum": kind, "signature": s, "exposure": e["w"], "ci_low": e["lo"], "ci_high": e["hi"],
                         "cosine": r["cosine"], "reads_accounted": r["reads_accounted"], "n_signatures": r["elbow"]})
pd.DataFrame(fit_rows, columns=["pool", "spectrum", "signature", "exposure", "ci_low", "ci_high", "cosine",
                                "reads_accounted", "n_signatures"]).to_csv("error_signatures.csv", index=False)

with open("versions.yml", "w") as fh:
    fh.write('"${task.process}":\\n')
    fh.write(f"    python: {platform.python_version()}\\n")
    fh.write(f"    pandas: {pd.__version__}\\n")
    fh.write(f"    numpy: {np.__version__}\\n")
