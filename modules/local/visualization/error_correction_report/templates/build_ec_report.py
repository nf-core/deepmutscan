#!/usr/bin/env python

# Build the self-contained error-correction HTML report for a whole run.
#
# One report covers every sequencing file: the page carries each file's numbers separately and
# recomputes every panel in the browser from whichever files are selected, so deselecting down to a
# single file reproduces exactly what a per-file report used to show.
#
# Nextflow template: values below are interpolated; backslashes meant for Python are doubled.

import base64
import json
import platform
import re

import numpy as np
import pandas as pd

sample = "$sample"
method = "$method"
# The false-doubles estimator actually used (mle|eb), so the report can name it. Empty/other for the
# wildtype method, where it does not apply.
estimator = "$false_doubles_method"
# Each entry is {id, type, replicate} - type/replicate come straight from the samplesheet (general
# for any protein, not hard-coded here) and drive the input-vs-output split in the report.
files = json.loads(base64.b64decode("$ids_b64").decode("utf-8"))
paths = "$corrected".split()

if len(files) != len(paths):
    raise SystemExit(f"{len(files)} ids but {len(paths)} files staged")

MUT = re.compile(r"^([A-Za-z*]+)(\\d+)([A-Za-z*]+)")

# `base_mut` holds the nucleotide change(s) of a variant as "<pos>:<ref>><alt>", comma-separated when
# the codon carries 2-3 changes (the NNK case). Positions are reference coordinates.
SNV = re.compile(r"^(\\d+):([ACGTacgt])>([ACGTacgt])\$")

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

# The six strand-symmetric substitution classes, pyrimidine-centric: a purine reference is
# reverse-complemented, so C>A and G>T are one and the same event. Order is canonical (as used for
# mutational signatures) and fixed - the colour of a class never depends on the data.
SNV_CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SNV_LABELS = {
    "C>A": "C>A / G>T", "C>G": "C>G / G>C", "C>T": "C>T / G>A",
    "T>A": "T>A / A>T", "T>C": "T>C / A>G", "T>G": "T>G / A>C",
}


def parse_pos(s):
    m = MUT.match(str(s))
    return (m.group(1), int(m.group(2)), m.group(3)) if m else (None, None, None)


def parse_snvs(base_mut):
    """"1220:T>G, 1221:T>G" -> [(1220,'T','G'), (1221,'T','G')]"""
    out = []
    for part in str(base_mut).split(","):
        m = SNV.match(part.strip())
        if m:
            out.append((int(m.group(1)), m.group(2).upper(), m.group(3).upper()))
    return out


def classify_snv(ref, alt):
    """Collapse a substitution onto its pyrimidine-centric class key."""
    if ref not in COMPLEMENT or alt not in COMPLEMENT or ref == alt:
        return None
    if ref in ("A", "G"):
        ref, alt = COMPLEMENT[ref], COMPLEMENT[alt]
    key = ref + ">" + alt
    return key if key in SNV_LABELS else None


def numcol(df, name):
    return (pd.to_numeric(df[name], errors="coerce") if name in df.columns
            else pd.Series([np.nan] * len(df)))


def f(x, nd):
    return None if x is None or x != x else round(float(x), nd)


# Variants are keyed on `codon_mut`, not `pos_mut`. `pos_mut` is NOT unique - one amino-acid
# substitution is reachable through up to three codons (6898 unique pos_mut across 10254 rows on the
# GID1A run), so keying on it would silently merge distinct nucleotide variants. `codon_mut` is
# unique (10254/10254).
#
# The files do not share a variant set either (union 10771 vs intersection 6886 on GID1A), so this is
# a union: a variant absent from a file gets None there, never 0. Absent means "not observed above
# threshold in that file", which is missing data, not a measured zero - averaging it in as 0 would
# drag every frequency down. The browser skips Nones and reports how many files contributed.
variants = {}
n_files = len(files)

for fi, (fid, path) in enumerate(zip(files, paths)):
    df = pd.read_csv(path)
    raw = numcol(df, "counts_raw")
    cor = numcol(df, "counts")
    cpc_raw = numcol(df, "counts_per_cov_raw")
    cpc_cor = numcol(df, "counts_per_cov")
    posmut = df["pos_mut"].astype(str)
    basemut = df["base_mut"].astype(str) if "base_mut" in df.columns else pd.Series([""] * len(df))
    codonmut = df["codon_mut"].astype(str) if "codon_mut" in df.columns else posmut

    for i in range(len(df)):
        r = raw.iloc[i]
        if pd.isna(r):
            continue
        key = codonmut.iloc[i].strip()
        v = variants.get(key)
        if v is None:
            wt, pos, _mut = parse_pos(posmut.iloc[i])
            nt = basemut.iloc[i].strip()
            snvs = parse_snvs(nt)
            # A class is only meaningful when a single nucleotide changed; NNK codon variants
            # (2-3 changes) mix classes and stay unclassified by design.
            kls = classify_snv(snvs[0][1], snvs[0][2]) if len(snvs) == 1 else None
            v = variants[key] = {
                "key": key, "pos_mut": posmut.iloc[i], "pos": pos, "wt_aa": wt,
                "nt": nt, "cls": kls,
                "ntpos": (snvs[0][0] if len(snvs) == 1 else None),
                "raw": [None] * n_files, "cor": [None] * n_files,
                "f_raw": [None] * n_files, "f_cor": [None] * n_files,
            }
        c = cor.iloc[i]
        v["raw"][fi] = f(r, 3)
        v["cor"][fi] = 0.0 if pd.isna(c) else f(c, 3)
        v["f_raw"][fi] = f(cpc_raw.iloc[i], 10)
        v["f_cor"][fi] = 0.0 if pd.isna(cpc_cor.iloc[i]) else f(cpc_cor.iloc[i], 10)

data = {
    "sample": sample,
    "method": method,
    "estimator": estimator,
    "files": files,
    "classes": [{"cls": k, "label": SNV_LABELS[k]} for k in SNV_CLASSES],
    "variants": sorted(variants.values(), key=lambda v: (v["pos"] if v["pos"] is not None else 0, v["key"])),
}

data_json = json.dumps(data, separators=(",", ":")).replace("</", "<\\\\/")

html = open("$template_html", encoding="utf-8").read()
html = html.replace("__EC_DATA__", data_json).replace("__SAMPLE__", sample)
with open("error_correction_report.html", "w", encoding="utf-8") as fh:
    fh.write(html)

with open("versions.yml", "w") as vf:
    vf.write('"${task.process}":\\n')
    vf.write("    python: " + platform.python_version() + "\\n")
    vf.write("    pandas: " + pd.__version__ + "\\n")
    vf.write("    numpy: " + np.__version__ + "\\n")

import os
mb = os.path.getsize("error_correction_report.html") / 1e6
n_snv = sum(1 for v in data["variants"] if v["cls"])
print(f"Wrote error_correction_report.html ({mb:.1f} MB) - {len(files)} files, "
      f"{len(data['variants'])} variants ({n_snv} single-nucleotide)")
