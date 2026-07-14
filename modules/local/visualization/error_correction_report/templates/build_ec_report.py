#!/usr/bin/env python

# Build the self-contained error-correction HTML report for one sample.
# Nextflow template: values below are interpolated; backslashes meant for Python are doubled.

import json
import re
import pandas as pd
import numpy as np

sample = "${meta.sample ?: meta.id}"
method = "$method"

df = pd.read_csv("$corrected")

raw = pd.to_numeric(df["counts_raw"], errors="coerce")
cor = pd.to_numeric(df["counts"], errors="coerce")
cov = pd.to_numeric(df["cov"], errors="coerce") if "cov" in df.columns else pd.Series([np.nan] * len(df))
posmut = df["pos_mut"].astype(str)

MUT = re.compile(r"^([A-Za-z*]+)(\\d+)([A-Za-z*]+)")

def parse_pos(s):
    m = MUT.match(s)
    if not m:
        return (None, None, None)
    return (m.group(1), int(m.group(2)), m.group(3))

variants = []
pos_agg = {}
pcts = []
for i in range(len(df)):
    s = posmut.iloc[i]
    wt, pos, _mut = parse_pos(s)
    r = raw.iloc[i]
    if pd.isna(r):
        continue
    r = float(r)
    c = 0.0 if pd.isna(cor.iloc[i]) else float(cor.iloc[i])
    delta = round(c - r, 3)
    pct = round((c - r) / r * 100.0, 2) if r > 0 else None
    cv = None if pd.isna(cov.iloc[i]) else int(cov.iloc[i])
    variants.append({"pos_mut": s, "pos": pos, "raw": round(r, 3), "corrected": round(c, 3),
                     "delta": delta, "pct": pct, "cov": cv})
    if pct is not None:
        pcts.append(pct)
    if pos is not None:
        a = pos_agg.setdefault(pos, {"pos": pos, "wt_aa": wt, "raw": 0.0, "cor": 0.0, "n": 0})
        a["raw"] += r
        a["cor"] += c
        a["n"] += 1

positions = []
for pos in sorted(pos_agg):
    a = pos_agg[pos]
    frac = (1.0 - a["cor"] / a["raw"]) if a["raw"] > 0 else None
    positions.append({"pos": pos, "wt_aa": a["wt_aa"],
                      "removed_frac": (round(frac, 6) if frac is not None else None),
                      "n": a["n"]})

# per-variant % change histogram (dynamic, symmetric-ish range)
hist = []
if pcts:
    arr = np.array(pcts, dtype=float)
    lo = float(np.floor(min(arr.min(), -1)))
    hi = float(np.ceil(max(arr.max(), 1)))
    counts, edges = np.histogram(arr, bins=18, range=(lo, hi))
    for k in range(len(counts)):
        hist.append({"x0": round(float(edges[k]), 3), "x1": round(float(edges[k + 1]), 3),
                     "count": int(counts[k])})

total_raw = float(np.nansum(raw.values))
total_cor = float(np.nansum(cor.values))
n_corrected = int(sum(1 for v in variants if v["delta"] != 0))
summary = {
    "n_variants": len(variants),
    "n_corrected": n_corrected,
    "total_raw": round(total_raw, 1),
    "total_corrected": round(total_cor, 1),
    "pct_removed": round((total_raw - total_cor) / total_raw * 100.0, 2) if total_raw > 0 else 0.0,
}

data = {"sample": sample, "method": method, "summary": summary,
        "positions": positions, "hist": hist, "variants": variants}

data_json = json.dumps(data, separators=(",", ":")).replace("</", "<\\\\/")

html = open("$template_html", encoding="utf-8").read()
html = html.replace("__EC_DATA__", data_json).replace("__SAMPLE__", sample)
with open("${meta.id}_error_correction_report.html", "w", encoding="utf-8") as fh:
    fh.write(html)

import platform
with open("versions.yml", "w") as vf:
    vf.write('"${task.process}":\\n')
    vf.write("    python: " + platform.python_version() + "\\n")
    vf.write("    pandas: " + pd.__version__ + "\\n")
    vf.write("    numpy: " + np.__version__ + "\\n")
