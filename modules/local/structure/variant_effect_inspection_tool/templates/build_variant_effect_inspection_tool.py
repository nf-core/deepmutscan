#!/usr/bin/env python3
"""
Build a single self-contained, interactive variant effect inspection tool (HTML) for nf-core/deepmutscan.

Nextflow ``template`` script for the VARIANT_EFFECT_INSPECTION_TOOL module (staged from this module's
``templates/`` folder). Reads the wildtype structure (PDB), the fitness estimation table and the
per-sample variant-count tables, aggregates per-residue/per-variant metrics, and writes one portable
HTML file with 3Dmol.js, the PDB coordinates and all data inlined (no external/CDN/sidecar
dependency), plus a structure_data.json sidecar and versions.yml.

Note: this is a Nextflow template, so input paths / the sample name are filled in as Nextflow
variables and every backslash meant for Python is written doubled (as in the pipeline's R templates).
"""
import json
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
MUT_RE = re.compile(r"^([A-Z])(\\d+)([A-Z*])\$")


def parse_pdb(pdb_path):
    """Return (pdb_text, residues) where residues = [{resnum, aa, plddt}] from CA atoms."""
    text = Path(pdb_path).read_text()
    residues = []
    seen = set()
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            resnum = int(line[22:26])
            bfac = float(line[60:66])
        except ValueError:
            continue
        if resnum in seen:
            continue
        seen.add(resnum)
        resname = line[17:20].strip().upper()
        residues.append({"resnum": resnum, "aa": THREE_TO_ONE.get(resname, "X"), "plddt": bfac})
    residues.sort(key=lambda r: r["resnum"])
    return text, residues


def best_alignment(wt_seq, struct_residues):
    """Find offset so that wt position i (1-based) maps to struct_residues[off + i - 1].

    Scans all offsets, maximising identity. Returns (offset0, identity_fraction) where the
    structure index for wt index j (0-based) is off0 + j.
    """
    struct_seq = "".join(r["aa"] for r in struct_residues)
    n, m = len(wt_seq), len(struct_seq)
    best_off, best_matches = 0, -1
    for off in range(-n + 1, m):
        matches = 0
        for j in range(n):
            k = off + j
            if 0 <= k < m and struct_seq[k] == wt_seq[j]:
                matches += 1
        if matches > best_matches:
            best_matches, best_off = matches, off
    identity = best_matches / max(1, n)
    return best_off, identity


def r6(x):
    """Round to 6 dp so output is byte-stable regardless of floating-point summation order."""
    return round(x, 6) if isinstance(x, (int, float)) else x


def to_float(x):
    try:
        f = float(x)
        return f if f == f else None  # drop NaN
    except (TypeError, ValueError):
        return None


def load_fitness(path):
    """Per-position variant records from the default fitness_estimation.tsv.

    Captures the mean fitness, its SD and the per-replicate (rescaled) fitness so the viewer can
    derive replicate-based metrics. Missing columns/values are tolerated (single-replicate data).
    """
    df = pd.read_csv(path, sep="\\t", dtype=str, keep_default_na=False)

    def col(name):
        return next((c for c in df.columns if c.strip().lower() == name), None)

    fit_col, sd_col = col("mean fitness"), col("fitness sd")
    r1_col, r2_col = col("rescaled_fitness_rep1"), col("rescaled_fitness_rep2")
    per_pos = {}
    for _, row in df.iterrows():
        pos = to_float(row.get("pos", ""))
        wt = (row.get("wt_aa", "") or "").strip()
        mut = (row.get("mut_aa", "") or "").strip()
        if pos is None or not wt or not mut:
            continue
        pos = int(pos)
        is_stop = mut == "*" or (row.get("stop", "").strip() not in ("", "0"))
        per_pos.setdefault(pos, {"wt": wt, "variants": []})["variants"].append({
            "mut": mut,
            "fitness": to_float(row.get(fit_col, "")) if fit_col else None,
            "sd": to_float(row.get(sd_col, "")) if sd_col else None,
            "rep1": to_float(row.get(r1_col, "")) if r1_col else None,
            "rep2": to_float(row.get(r2_col, "")) if r2_col else None,
            "n_out": (to_float(row.get("output1", "")) or 0) + (to_float(row.get("output2", "")) or 0),
            "syn": mut == wt, "stop": bool(is_stop),
        })
    return per_pos


def load_counts(paths):
    """Aggregate coverage / counts / counts_per_cov across all sample CSVs.

    Returns two dicts: per position ``{pos: {...}}`` (used to colour the structure) and per variant
    ``{(pos, mut): {...}}`` (used for the substitution-specific bee-swarm distributions).
    """
    frames = []
    for p in paths:
        try:
            frames.append(pd.read_csv(p))
        except Exception as exc:  # pragma: no cover - defensive
            print(f"WARN: could not read {p}: {exc}", file=sys.stderr)
    if not frames:
        return {}, {}
    df = pd.concat(frames, ignore_index=True)
    # The "A1M"-style WT/pos/mut label lives in the `pos_mut` column
    # (`aa_mut` holds e.g. "M:A>M"). Fall back to aa_mut only if pos_mut is absent.
    label_col = "pos_mut" if "pos_mut" in df.columns else ("aa_mut" if "aa_mut" in df.columns else None)
    if label_col is None:
        return {}, {}
    cols = (("coverage", "cov"), ("counts", "counts"), ("counts_per_cov", "counts_per_cov"),
            ("counts_raw", "counts_raw"))
    pos_acc, var_acc = {}, {}
    for _, row in df.iterrows():
        m = MUT_RE.match(str(row.get(label_col, "")).strip())
        if not m:
            continue
        pos, mut = int(m.group(2)), m.group(3)
        pa = pos_acc.setdefault(pos, {k: [] for k, _ in cols})
        va = var_acc.setdefault((pos, mut), {k: [] for k, _ in cols})
        for key, src in cols:
            v = to_float(row.get(src, ""))
            if v is not None:
                pa[key].append(v)
                va[key].append(v)

    def meaned(acc):
        return {k: {kk: (sum(vv) / len(vv) if vv else None) for kk, vv in a.items()}
                for k, a in acc.items()}

    return meaned(pos_acc), meaned(var_acc)


class args:
    # Filled in by Nextflow when this template is staged and run by the VARIANT_EFFECT_INSPECTION_TOOL process.
    pdb = "$pdb"
    fitness = "$fitness_tsv"
    aa_seq = "$aa_seq"
    counts = "$variant_counts".split()          # collected per-sample CSVs (staged as vc*/*)
    threedmol = "$threedmol_js"
    template = "$template_html"
    sample = "${meta.sample ?: meta.id}"
    out_html = "variant_effect_inspection_tool.html"
    out_json = "structure_data.json"


def main():

    wt_seq = Path(args.aa_seq).read_text().strip().replace("\\n", "").replace("\\r", "")
    pdb_text, struct_residues = parse_pdb(args.pdb)
    off0, identity = best_alignment(wt_seq, struct_residues)
    struct_by_index = struct_residues  # 0-based list

    fit_per_pos = load_fitness(args.fitness)
    # per-position aggregate (structure colour) and per-variant values (substitution-specific swarm)
    counts_per_pos, counts_per_var = load_counts(sorted(args.counts))  # sorted -> deterministic

    # ---- candidate metrics (dropped automatically when a dataset has no data for them) ----
    # scope="variant": has a per-substitution value (used in the bee-swarm); the position mean is
    # still what colours the structure. reverse=True flips the colour ramp (low -> red, high -> blue).
    METRIC_SPECS = [
        ("fitness", dict(label="Fitness (nf-core/deepmutscan)", type="diverging", scope="variant",
                         reverse=True,
                         desc="Fitness per substitution; the position mean is shown on the structure")),
        ("fitness_sd", dict(label="Fitness uncertainty (SD)", type="sequential", scope="variant",
                            desc="Standard deviation of the fitness estimate (lower = more confident)")),
        ("accuracy", dict(label="Prediction accuracy (pLDDT)", type="plddt", scope="residue", fixed=[0, 100],
                          desc="Per-residue model confidence from the structure B-factor")),
        ("coverage", dict(label="Coverage", type="sequential", scope="variant",
                          desc="Sequencing coverage of this residue (the same for all substitutions at one residue)")),
        ("counts", dict(label="Counts", type="sequential", scope="variant",
                        desc="Variant read counts per substitution (position mean shown on the structure)")),
        ("counts_per_cov", dict(label="Counts / coverage", type="sequential", scope="variant",
                                desc="Counts normalised by coverage, per substitution (position mean shown on the structure)")),
        ("n_substitutions", dict(label="Measured substitutions", type="sequential", scope="residue",
                                 desc="Number of substitutions with a measured fitness at the position")),
        ("error_bias", dict(label="Positional error bias", type="sequential", scope="residue",
                            desc="Fraction of counts removed by sequencing-error correction at this position")),
    ]
    variant_metrics = {k for k, s in METRIC_SPECS if s["scope"] == "variant"}

    def mean(xs):
        return (sum(xs) / len(xs)) if xs else None

    def rep_corr(vs):
        pairs = [(v["rep1"], v["rep2"]) for v in vs
                 if v.get("rep1") is not None and v.get("rep2") is not None]
        if len(pairs) < 3:
            return None
        a = np.array([p[0] for p in pairs], float)
        b = np.array([p[1] for p in pairs], float)
        if a.std() == 0 or b.std() == 0:
            return None
        return float(np.corrcoef(a, b)[0, 1])

    residues = []
    for i, wt_aa in enumerate(wt_seq):  # i is 0-based; position = i + 1
        pos = i + 1
        k = off0 + i
        srec = struct_by_index[k] if 0 <= k < len(struct_by_index) else None
        variants = fit_per_pos.get(pos, {}).get("variants", [])
        miss = [v for v in variants if not v["syn"] and not v["stop"]]
        fit_vals = [v["fitness"] for v in miss if v["fitness"] is not None]
        sd_vals = [v["sd"] for v in miss if v["sd"] is not None]
        cp = counts_per_pos.get(pos, {})
        craw, cnow = cp.get("counts_raw"), cp.get("counts")
        error_bias = (1.0 - cnow / craw) if (craw not in (None, 0) and cnow is not None) else None
        res_metrics = {
            "fitness": r6(mean(fit_vals)),
            "fitness_sd": r6(mean(sd_vals)),
            "accuracy": r6(srec["plddt"]) if srec else None,
            "coverage": r6(cp.get("coverage")),
            "counts": r6(cp.get("counts")),
            "counts_per_cov": r6(cp.get("counts_per_cov")),
            "n_substitutions": (len(fit_vals) or None),
            "error_bias": r6(error_bias),
        }
        vlist = []
        for v in sorted(variants, key=lambda v: AA_ORDER.index(v["mut"]) if v["mut"] in AA_ORDER else 99):
            vc = counts_per_var.get((pos, v["mut"]), {})
            vlist.append({
                "mut": v["mut"], "syn": v["syn"], "stop": v["stop"], "n_out": r6(v["n_out"]),
                "metrics": {
                    "fitness": r6(v["fitness"]), "fitness_sd": r6(v["sd"]),
                    "coverage": r6(vc.get("coverage")), "counts": r6(vc.get("counts")),
                    "counts_per_cov": r6(vc.get("counts_per_cov")),
                },
            })
        residues.append({
            "pos": pos, "wt": wt_aa,
            "resnum": srec["resnum"] if srec else None,
            "mapped": bool(srec and srec["aa"] == wt_aa),
            "metrics": res_metrics, "variants": vlist,
        })

    # keep only metrics that actually have data (generality across diverse datasets)
    metrics = {}
    for key, spec in METRIC_SPECS:
        vals = [r["metrics"].get(key) for r in residues if r["metrics"].get(key) is not None]
        if not vals:
            continue
        if "fixed" in spec:
            mrange = spec["fixed"]
        elif spec["type"] == "diverging":
            a = max(abs(min(vals)), abs(max(vals))) or 1
            mrange = [r6(-a), r6(a)]
        else:
            lo, hi = min(vals), max(vals)
            mrange = [r6(lo), r6(hi if hi != lo else lo + 1)]
        metrics[key] = {"label": spec["label"], "type": spec["type"], "scope": spec["scope"],
                        "range": mrange, "desc": spec["desc"]}
        if spec.get("reverse"):
            metrics[key]["reverse"] = True

    present = set(metrics)
    for r in residues:
        r["metrics"] = {k: v for k, v in r["metrics"].items() if k in present}
        for v in r["variants"]:
            v["metrics"] = {k: val for k, val in v["metrics"].items()
                            if k in present and k in variant_metrics}

    data = {
        "sample": args.sample,
        "sequence": wt_seq,
        "alignment": {"offset": off0, "identity": round(identity, 4),
                      "n_residues_structure": len(struct_residues)},
        "metrics": metrics,
        "residues": residues,
    }

    Path(args.out_json).write_text(json.dumps(data, indent=2))

    threedmol_js = Path(args.threedmol).read_text()
    template = Path(args.template).read_text()
    html = (template
            .replace("/*__THREEDMOL_JS__*/", threedmol_js)
            .replace('"__PDB_TEXT__"', json.dumps(pdb_text))
            .replace("/*__DATA_JSON__*/", json.dumps(data))
            .replace("__SAMPLE__", args.sample))
    Path(args.out_html).write_text(html)

    Path("versions.yml").write_text(
        '"${task.process}":\\n'
        + "    python: " + platform.python_version() + "\\n"
        + "    pandas: " + pd.__version__ + "\\n"
        + "    numpy: " + np.__version__ + "\\n"
    )
    print(f"Wrote {args.out_html} (alignment offset={off0}, identity={identity:.2%}, "
          f"{len(residues)} residues)")


if __name__ == "__main__":
    main()
