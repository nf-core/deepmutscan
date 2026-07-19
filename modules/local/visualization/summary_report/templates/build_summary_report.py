#!/usr/bin/env python3
"""Build the self-contained all-in-one HTML report for nf-core/deepmutscan.

Nextflow ``template`` script for the VISUALIZATION_SUMMARY_REPORT module. Everything the report
shows arrives through one generic mechanism: a base64 JSON manifest describing each staged asset
(`kind|group|name`) plus the assets themselves, in the same order. Sections that have no entries
are simply omitted, so the same code serves runs with or without fitness / error correction / a PDB.

Two kinds of payload end up in the HTML:
  * data  - count tables and the fitness table, aggregated down to what the plots actually need,
            then rendered natively as interactive SVG in the browser.
  * blobs - the original R-generated plot PDFs and the standalone detail reports, embedded as
            base64 data URIs so the single file needs no sidecars to view, download or share.

Note: this is a Nextflow template, so every backslash meant for Python is doubled and literal
dollar signs are escaped.
"""
import base64
import json
import platform
import re
from pathlib import Path

import numpy as np
import pandas as pd

MUT = re.compile(r"^([A-Za-z*]+)(\\d+)([A-Za-z*]+)")
REP_RE = re.compile(r"^rescaled_fitness_rep(\\d+)\$", re.I)
IN_RE = re.compile(r"^input(\\d+)\$", re.I)
OUT_RE = re.compile(r"^output(\\d+)\$", re.I)

# Reference data for the citation page. Every entry is transcribed from CITATIONS.md, which is the
# pipeline's own maintained list - the point is that the report and the repo can never drift apart,
# and that nothing here is invented. `pkgs` are the names these tools use in the run's versions.yml,
# which is what decides whether a tool actually ran and therefore appears.
CITATIONS = [
    {"id": "deepmutscan", "always": True, "pkgs": [],
     "authors": "Wehnert B, et al.", "title": "nf-core/deepmutscan",
     "journal": "", "year": "", "key": "deepmutscan",
     "note": "bioRxiv preprint, coming soon."},
    {"id": "nfcore", "always": True, "pkgs": [],
     "authors": "Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S",
     "title": "The nf-core framework for community-curated bioinformatics pipelines",
     "journal": "Nature Biotechnology", "year": "2020", "volume": "38", "number": "3", "pages": "276--278",
     "doi": "10.1038/s41587-020-0439-x", "key": "Ewels2020"},
    {"id": "nextflow", "always": True, "pkgs": ["nextflow"],
     "authors": "Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C",
     "title": "Nextflow enables reproducible computational workflows",
     "journal": "Nature Biotechnology", "year": "2017", "volume": "35", "number": "4", "pages": "316--319",
     "doi": "10.1038/nbt.3820", "key": "DiTommaso2017"},
    {"id": "fastqc", "pkgs": ["fastqc"],
     "authors": "Andrews S", "title": "FastQC: A Quality Control Tool for High Throughput Sequence Data",
     "year": "2010", "url": "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/", "key": "Andrews2010"},
    {"id": "multiqc", "pkgs": ["multiqc"],
     "authors": "Ewels P, Magnusson M, Lundin S, K\\u00e4ller M",
     "title": "MultiQC: summarize analysis results for multiple tools and samples in a single report",
     "journal": "Bioinformatics", "year": "2016", "volume": "32", "number": "19", "pages": "3047--3048",
     "doi": "10.1093/bioinformatics/btw354", "key": "Ewels2016"},
    {"id": "samtools", "pkgs": ["samtools"],
     "authors": "Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H",
     "title": "Twelve years of SAMtools and BCFtools",
     "journal": "GigaScience", "year": "2021", "volume": "10", "number": "2", "pages": "giab008",
     "doi": "10.1093/gigascience/giab008", "key": "Danecek2021"},
    # BWA-MEM, not BWA-backtrack: the pipeline runs `bwa mem`, whose citation is the 2013 arXiv
    # preprint - the same one nf-core's own bwa/mem module points at. Li & Durbin 2009 describes
    # `bwa aln`, a different algorithm this pipeline never calls.
    {"id": "bwa", "pkgs": ["bwa"],
     "authors": "Li H", "title": "Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM",
     "journal": "arXiv", "year": "2013", "eprint": "1303.3997", "primaryclass": "q-bio.GN",
     "url": "https://arxiv.org/abs/1303.3997", "key": "Li2013"},
    {"id": "pysam", "pkgs": ["pysam"],
     "authors": "Heger A, et al.", "title": "pysam: a Python module for reading, manipulating and writing genomic data sets",
     "year": "", "url": "https://github.com/pysam-developers/pysam", "key": "pysam"},
    {"id": "polars", "pkgs": ["polars"],
     "authors": "Vink R, et al.", "title": "Polars: Lightning-fast DataFrame library for Rust and Python",
     "year": "", "url": "https://github.com/pola-rs/polars", "key": "polars"},
    {"id": "biopython", "pkgs": ["biopython"],
     "authors": "Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL",
     "title": "Biopython: freely available Python tools for computational molecular biology and bioinformatics",
     "journal": "Bioinformatics", "year": "2009", "volume": "25", "number": "11", "pages": "1422--1423",
     "doi": "10.1093/bioinformatics/btp163", "key": "Cock2009"},
    {"id": "threedmol", "pkgs": ["3dmol"],
     "authors": "Rego N, Koes D", "title": "3Dmol.js: molecular visualization with WebGL",
     "journal": "Bioinformatics", "year": "2015", "volume": "31", "number": "8", "pages": "1322--1324",
     "doi": "10.1093/bioinformatics/btu829", "key": "Rego2015",
     "note": "Bundled in the interactive variant effect inspection tool (BSD-3-Clause)."},
    {"id": "dimsum", "pkgs": ["r-dimsum", "dimsum"],
     "authors": "Faure AJ, Schmiedel JM, Baeza-Centurion P, Lehner B",
     "title": "DiMSum: an error model and pipeline for analyzing deep mutational scanning data and diagnosing common experimental pathologies",
     "journal": "Genome Biology", "year": "2020", "volume": "21", "number": "1", "pages": "207",
     "doi": "10.1186/s13059-020-02091-3", "key": "Faure2020"},
    {"id": "mutscan", "pkgs": ["r-mutscan", "mutscan"],
     "authors": "Soneson C, Bendel AM, Diss G, Stadler MB",
     "title": "mutscan - a flexible R package for efficient end-to-end analysis of multiplexed assays of variant effect data",
     "journal": "Genome Biology", "year": "2023", "volume": "24", "number": "1", "pages": "132",
     "doi": "10.1186/s13059-023-02967-0", "key": "Soneson2023"},
]


class args:
    run = "$run_b64"
    manifest = "$manifest_b64"
    assets = "$assets".split()
    template = "$template_html"
    threedmol = "$threedmol_js"          # 3Dmol.js, inlined only when a structure is shown
    out_html = "deepmutscan_report.html"


def b64json(s):
    return json.loads(base64.b64decode(s).decode("utf-8"))


def data_uri(path, mime):
    return "data:" + mime + ";base64," + base64.b64encode(Path(path).read_bytes()).decode("ascii")


# Embedded reports run in a frame whose origin is opaque, and in an opaque origin `localStorage` and
# `document.cookie` do not return empty - they *throw* SecurityError. MultiQC v1.35 reads both while
# bootstrapping without a try/catch, so the exception unwinds its ready handler and the report is left
# frozen on its "Loading report.." placeholder with no plots. Handing it a memory-backed store that
# never throws is enough; the only thing lost is persistence of the reader's MultiQC preferences.
# This prepends to the file we embed and does not touch the published multiqc_report.html.
STORAGE_SHIM = """<script>
(function () {
  try { window.localStorage.getItem('__probe__'); } catch (e) {
    var m = Object.create(null);
    var s = {
      getItem: function (k) { return k in m ? m[k] : null; },
      setItem: function (k, v) { m[k] = String(v); },
      removeItem: function (k) { delete m[k]; },
      clear: function () { m = Object.create(null); },
      key: function (i) { return Object.keys(m)[i] || null; }
    };
    Object.defineProperty(s, 'length', { get: function () { return Object.keys(m).length; } });
    try { Object.defineProperty(window, 'localStorage', { value: s, configurable: true }); } catch (e2) {}
  }
  try { void document.cookie; } catch (e) {
    var jar = '';
    try {
      Object.defineProperty(document, 'cookie', {
        get: function () { return jar; },
        set: function (v) { jar = String(v).split(';')[0]; },
        configurable: true
      });
    } catch (e2) {}
  }
})();
</script>"""

HEAD_RE = re.compile(r"<head[^>]*>", re.I)


def html_for_frame(path):
    """Read an embedded report and make it survive an opaque origin."""
    text = Path(path).read_text(encoding="utf-8", errors="surrogateescape")
    m = HEAD_RE.search(text)
    if m:
        text = text[: m.end()] + STORAGE_SHIM + text[m.end():]
    return text


def r6(x):
    return round(float(x), 6) if x is not None and x == x else None


def to_float(x):
    try:
        v = float(x)
        return round(v, 6) if v == v else None
    except (TypeError, ValueError):
        return None


def wt_orf(fasta_path, reading_frame):
    """Wildtype ORF nucleotides plus the reference offset they start at.

    --reading_frame is the 1-based inclusive range the pipeline uses everywhere else (e.g.
    "352-1383"), converted to a 0-based half-open slice here. `start` is kept so nucleotide positions
    can be reported in reference coordinates, which is what base_mut uses.

    Returns None rather than guessing whenever anything fails to line up: a wrong offset would
    mislabel every nucleotide change while still looking entirely plausible.
    """
    if not fasta_path or not reading_frame:
        return None
    try:
        lo_s, hi_s = str(reading_frame).split("-")
        lo, hi = int(lo_s), int(hi_s)
    except (TypeError, ValueError):
        return None
    seq = "".join(l.strip() for l in Path(fasta_path).read_text().splitlines()
                  if l.strip() and not l.startswith(">")).upper()
    if lo < 1 or hi > len(seq) or hi < lo:
        return None
    orf = seq[lo - 1:hi]
    if len(orf) % 3:
        return None
    return {"seq": orf, "start": lo}


def parse_pos(s):
    m = MUT.match(str(s))
    return (m.group(1), int(m.group(2)), m.group(3)) if m else (None, None, None)


def num(df, name):
    return (pd.to_numeric(df[name], errors="coerce") if name in df.columns
            else pd.Series([np.nan] * len(df)))


def per_position(path, orf_len):
    """Aggregate a filtered-by-library table to per-ORF-position series.

    Deliberately mirrors `global_position_biases_{counts,cov}.R` so the interactive curves and the
    pipeline's published PDFs describe the same quantity:

      * coverage  - a property of the position, so it is averaged over the substitutions measured
                    there (the R comment claims it is constant; on real data it varies by ~0.25%
                    because reads must span every varying base plus --min_flank, which is small
                    enough that the mean is the right summary).
      * counts / counts_per_cov
                  - averaged *within* a varying-bases class and never pooled across classes. The
                    classes differ several-fold in both counts and frequency, and their mix changes
                    along the ORF, so a pooled curve would track composition rather than the library.
      * R uses na.rm=FALSE for these per-position means, hence skipna=False here.

    Series run over every ORF position 1..orf_len, with None where a class has no variants, so that
    the rolling mean in the browser sees the same gaps rollapply() does.

    `counts_raw` is only present once error correction has run, and drives the per-sample summary.
    """
    df = pd.read_csv(path)
    pos = df["pos_mut"].astype(str).map(lambda s: parse_pos(s)[1])
    cov, cnt = num(df, "cov"), num(df, "counts")
    cpc, raw = num(df, "counts_per_cov"), num(df, "counts_raw")
    vb = num(df, "varying_bases")
    d = pd.DataFrame({"pos": pos, "cov": cov, "counts": cnt, "cpc": cpc, "raw": raw, "vb": vb})
    d = d.dropna(subset=["pos"])
    d["pos"] = d["pos"].astype(int)

    n = int(orf_len) if orf_len else (int(d["pos"].max()) if len(d) else 0)
    axis = list(range(1, n + 1))

    cov_by = d.groupby("pos")["cov"].mean(numeric_only=False)
    positions = [{"pos": p, "cov": r6(cov_by.get(p))} for p in axis]

    classes = {}
    for k in (1, 2, 3):
        sub = d[d["vb"] == k]
        cm = sub.groupby("pos")["counts"].apply(lambda s: s.mean(skipna=False))
        fm = sub.groupby("pos")["cpc"].apply(lambda s: s.mean(skipna=False))
        classes[str(k)] = [{"pos": p, "counts": r6(cm.get(p)), "cpc": r6(fm.get(p))} for p in axis]

    ec = None
    if raw.notna().any():
        tr, tc = float(np.nansum(raw.values)), float(np.nansum(cnt.values))
        ec = {
            "n_variants": int(len(d)),
            "n_corrected": int((raw.fillna(0) != cnt.fillna(0)).sum()),
            "total_raw": round(tr, 1),
            "total_corrected": round(tc, 1),
            "pct_removed": round((tr - tc) / tr * 100.0, 2) if tr > 0 else 0.0,
        }
    return positions, classes, ec


def heatmap(path):
    """Counts-per-coverage grid (position x substituted amino acid) for one sample."""
    df = pd.read_csv(path)
    if "position" not in df.columns or "mut_aa" not in df.columns:
        return None
    pos = pd.to_numeric(df["position"], errors="coerce")
    val = num(df, "total_counts_per_cov")
    cells = [{"p": int(p), "a": str(a), "v": r6(v)}
             for p, a, v in zip(pos, df["mut_aa"].astype(str), val)
             if p == p and v == v]
    wt = {}
    if "wt_aa" in df.columns:
        for p, w in zip(pos, df["wt_aa"].astype(str)):
            if p == p and w and w != "nan":
                wt[int(p)] = w
    return {"cells": cells, "wt": wt}


def codon_changes(nt_seq_field, pos, wt_nt):
    """Nucleotide changes behind one amino-acid substitution.

    `nt_seq` holds the full ORF for every distinct coding variant that produces this amino-acid
    change, comma-separated - an NNK library typically reaches one residue through several codons.
    Each is diffed against the wildtype ORF over this residue's codon only: changes elsewhere belong
    to other residues and would be double counted. Returns one entry per distinct codon.

    Reported positions are in REFERENCE coordinates, not ORF-relative ones, because that is the
    convention `base_mut` uses everywhere else in the pipeline - checked against variantCounts, where
    4000/4000 variants agree with the reference convention and 0/4000 with the ORF-relative one.
    Emitting ORF-relative numbers here would be off by the reading frame's start and would silently
    disagree with the error-correction report.
    """
    if not wt_nt or not wt_nt.get("seq") or not nt_seq_field:
        return []
    seq_wt, rf_lo = wt_nt["seq"], wt_nt["start"]
    lo = (int(pos) - 1) * 3
    if lo + 3 > len(seq_wt):
        return []
    wt_codon = seq_wt[lo:lo + 3]
    out, seen = [], set()
    for seq in (s.strip() for s in str(nt_seq_field).split(",")):
        if len(seq) < lo + 3:
            continue
        mut_codon = seq[lo:lo + 3]
        if mut_codon == wt_codon or mut_codon in seen:
            continue
        seen.add(mut_codon)
        changes, positions = [], []
        for k in range(3):
            if wt_codon[k] != mut_codon[k]:
                changes.append(f"{wt_codon[k]}>{mut_codon[k]}")
                positions.append(rf_lo + lo + k)   # reference coords, as base_mut reports them
        out.append({"wt_codon": wt_codon, "mut_codon": mut_codon,
                    "changes": changes, "positions": positions})
    return out


def load_fitness(path, wt_nt=None):
    """Per-substitution fitness plus its per-position mean, from the default fitness table."""
    df = pd.read_csv(path, sep="\\t", dtype=str, keep_default_na=False)

    def col(n):
        return next((c for c in df.columns if c.strip().lower() == n), None)

    fc, sc = col("mean fitness"), col("fitness sd")
    if not fc:
        return None
    # Replicate count is a per-run property, so discover the columns instead of assuming rep1/rep2.
    rep_cols = sorted((c for c in df.columns if REP_RE.match(c.strip())),
                      key=lambda c: int(REP_RE.match(c.strip()).group(1)))
    in_cols = [c for c in df.columns if IN_RE.match(c.strip())]
    out_cols = [c for c in df.columns if OUT_RE.match(c.strip())]
    ntc = col("nt_seq")
    rows, per_pos = [], {}
    for _, r in df.iterrows():
        try:
            p = int(float(r.get("pos", "") or "nan"))
        except (TypeError, ValueError):
            continue
        wt, mut = (r.get("wt_aa", "") or "").strip(), (r.get("mut_aa", "") or "").strip()
        if not wt or not mut:
            continue
        try:
            f = float(r.get(fc, ""))
        except (TypeError, ValueError):
            continue
        if f != f:
            continue
        sd = None
        try:
            sd = float(r.get(sc, "")) if sc else None
            if sd != sd:
                sd = None
        except (TypeError, ValueError):
            sd = None
        stop = mut == "*" or (r.get("stop", "").strip() not in ("", "0"))

        def tot(cols):
            s = 0.0
            got = False
            for c in cols:
                try:
                    s += float(r.get(c, "") or 0)
                    got = True
                except (TypeError, ValueError):
                    pass
            return s if got else None

        cell = {"p": p, "wt": wt, "a": mut, "f": r6(f), "sd": r6(sd),
                "syn": mut == wt, "stop": bool(stop)}
        if rep_cols:
            cell["reps"] = [to_float(r.get(c, "")) for c in rep_cols]
        cell["n_in"], cell["n_out"] = tot(in_cols), tot(out_cols)
        if ntc and wt_nt:
            nt = codon_changes(r.get(ntc, ""), p, wt_nt)
            if nt:
                cell["nt"] = nt
        rows.append(cell)
        if mut != wt and not stop:
            per_pos.setdefault(p, []).append(f)
    return {"cells": rows}


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_pdb(pdb_path):
    """(pdb_text, residues) with residues = [{resnum, aa}] from the CA atoms, ordered by residue number."""
    text = Path(pdb_path).read_text()
    residues, seen = [], set()
    for line in text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        if resnum in seen:
            continue
        seen.add(resnum)
        residues.append({"resnum": resnum, "aa": THREE_TO_ONE.get(line[17:20].strip().upper(), "X")})
    residues.sort(key=lambda r: r["resnum"])
    return text, residues


def best_alignment(wt_seq, struct_residues):
    """Offset mapping wt position j (0-based) to struct_residues[off0 + j], maximising identity.

    The PDB's residue numbering need not start at the ORF's first residue (AFDB models often carry a
    signal peptide or a different frame), so the mapping is found by sliding the sequences, exactly as
    the inspection tool does - the two must agree or the structure would be coloured off by a frame.
    """
    struct_seq = "".join(r["aa"] for r in struct_residues)
    n, m = len(wt_seq), len(struct_seq)
    best_off, best_matches = 0, -1
    for off in range(-n + 1, m):
        matches = sum(1 for j in range(n) if 0 <= off + j < m and struct_seq[off + j] == wt_seq[j])
        if matches > best_matches:
            best_matches, best_off = matches, off
    return best_off, best_matches / max(1, n)


def structure_data(pdb_path, wt_seq, fitness):
    """PDB text plus per-residue mean fitness, keyed by PDB residue number, for the spinning viewer.

    The colour at a residue is the mean fitness over the missense substitutions measured there - the
    same position summary the fitness heatmap collapses to - so the structure and the heatmap tell one
    story. Returns None (viewer simply hidden) whenever a structure cannot be trusted: no PDB, no
    fitness, or an alignment that never lines up.
    """
    if not pdb_path or not wt_seq or not fitness or not fitness.get("cells"):
        return None
    pdb_text, residues = parse_pdb(pdb_path)
    if not residues:
        return None
    off0, identity = best_alignment(wt_seq, residues)
    by_pos = {}
    for c in fitness["cells"]:
        if c.get("syn") or c.get("stop") or c.get("f") is None:
            continue
        by_pos.setdefault(c["p"], []).append(c["f"])
    pos_mean = {p: sum(v) / len(v) for p, v in by_pos.items() if v}
    resi_fit = {}
    for j in range(len(wt_seq)):
        k = off0 + j
        if 0 <= k < len(residues) and (j + 1) in pos_mean:
            resi_fit[residues[k]["resnum"]] = round(pos_mean[j + 1], 4)
    if not resi_fit:
        return None
    fmax = max(abs(v) for v in resi_fit.values()) or 1.0
    return {"resi_fit": resi_fit, "fmax": round(fmax, 4),
            "identity": round(identity, 4), "n_residues": len(residues), "pdb": pdb_text}


def _dur_s(s):
    """Nextflow-formatted duration ('7m 37s', '746ms', '1h 2m') to seconds."""
    s = str(s).strip()
    if not s or s == "-":
        return 0.0
    tot = 0.0
    for num, unit in re.findall(r"([\\d.]+)\\s*(ms|s|m|h|d)", s):
        tot += float(num) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}[unit]
    return tot


def _size_b(s):
    """Nextflow-formatted size ('3 GB', '4.9 MB') to bytes."""
    m = re.match(r"([\\d.]+)\\s*(B|KB|MB|GB|TB)", str(s).strip(), re.I)
    if not m:
        return 0.0
    return float(m.group(1)) * {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}[m.group(2).upper()]


def run_stats(path):
    """Per-process compute stats from Nextflow's execution trace.

    The trace is written continuously, so at report time it holds every task finished so far - i.e.
    everything but the report itself and the version collation that follows it. That is the whole
    compute-heavy body of the run, which is what the page is about; it is labelled 'as of report
    generation' rather than pretended to be the final total. CACHED tasks are counted separately
    because they consumed no compute this run.
    """
    df = pd.read_csv(path, sep="\\t", dtype=str, keep_default_na=False)
    if "name" not in df.columns:
        return None

    def proc(n):
        return re.sub(r"\\s*\\(.*\\)\\s*\$", "", n).split(":")[-1]  # drop "(sample)" and the module path

    rows = {}
    for _, r in df.iterrows():
        p = proc(r.get("name", ""))
        d = rows.setdefault(p, {"n": 0, "cached": 0, "cpu_time": 0.0, "rss": 0.0, "cpu": []})
        d["n"] += 1
        if r.get("status", "") == "CACHED":
            d["cached"] += 1
        d["cpu_time"] += _dur_s(r.get("realtime", ""))
        d["rss"] = max(d["rss"], _size_b(r.get("peak_rss", "")))
        c = r.get("%cpu", "").replace("%", "").strip()
        try:
            d["cpu"].append(float(c))
        except ValueError:
            pass
    procs = [{
        "process": p, "tasks": d["n"], "cached": d["cached"],
        "cpu_time": round(d["cpu_time"], 1), "peak_rss": round(d["rss"]),
        "cpu_pct": round(sum(d["cpu"]) / len(d["cpu"])) if d["cpu"] else None,
    } for p, d in rows.items()]
    procs.sort(key=lambda x: -x["cpu_time"])
    status = df["status"] if "status" in df.columns else pd.Series([], dtype=str)
    total = {
        "tasks": int(df.shape[0]),
        "cached": int((status == "CACHED").sum()),
        "cpu_time": round(sum(p["cpu_time"] for p in procs), 1),
        "peak_rss": round(max((p["peak_rss"] for p in procs), default=0)),
    }
    return {"processes": procs, "total": total}


def software_versions(path):
    """Parse pipeline_info's software versions YAML into {process: {package: version}}.

    Hand-parsed rather than via a YAML library: the file is two levels of trivial `key: value` and
    pandas' container carries no yaml. Duplicate process keys are merged rather than overwritten -
    they occur in practice when a module fails to interpolate its own name.
    """
    out = {}
    cur = None
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if not line.startswith((" ", "\\t")):
            cur = line.strip().rstrip(":").strip('"')
            out.setdefault(cur, {})
        elif cur:
            k, _, v = line.strip().partition(":")
            if k:
                out[cur][k.strip()] = v.strip()
    return out


def citations(versions, extra=()):
    """The subset of CITATIONS whose tools actually ran, plus the pipeline and frameworks always.

    Package names in versions.yml are the primary signal, but they are not always sufficient: DiMSum,
    for instance, reports only `r-base`, so it is invisible at the package level. `extra` carries the
    tokens that a run's own flags settle unambiguously (--dimsum, --mutscan, a --pdb structure). A tool
    outside both signals is never cited - the page is about this run, not the pipeline's full
    repertoire.
    """
    seen = {p.strip().lower() for pkgs in versions.values() for p in pkgs}
    seen |= {e.lower() for e in extra}
    out = []
    for c in CITATIONS:
        if c.get("always") or any(p.lower() in seen for p in c["pkgs"]):
            out.append({k: v for k, v in c.items() if k != "pkgs"})
    return out


def biblatex(entries):
    """A biblatex source block for the given entries."""
    lines = []
    for c in entries:
        kind = "@article"
        if c.get("eprint"):
            kind = "@misc"
        elif not c.get("journal"):
            kind = "@misc"
        lines.append(f"{kind}{{{c['key']},")
        lines.append("  author       = {" + c["authors"].replace(", ", " and ") + "},")
        lines.append("  title        = {{" + c["title"] + "}},")
        for fld, key in (("journaltitle", "journal"), ("year", "year"), ("volume", "volume"),
                         ("number", "number"), ("pages", "pages"), ("doi", "doi"),
                         ("eprint", "eprint"), ("eprintclass", "primaryclass"), ("url", "url")):
            if c.get(key):
                lines.append(f"  {fld:<12} = {{{c[key]}}},")
        if c.get("note"):
            lines.append("  note         = {" + c["note"] + "},")
        lines.append("}")
        lines.append("")
    return "\\n".join(lines)


def main():
    run = b64json(args.run)
    manifest = b64json(args.manifest)
    if len(manifest) != len(args.assets):
        raise SystemExit(f"manifest has {len(manifest)} entries but {len(args.assets)} files were staged")

    samples, plots, reports, fitness, dl, versions = {}, {}, [], None, [], {}

    def sample(sid):
        return samples.setdefault(sid, {"id": sid, "positions": [], "classes": {}, "heatmap": None, "ec": None, "seqdepth": None})

    # The ORF length has to come from the wildtype protein sequence, not from the highest position that
    # happens to carry a variant, or trailing residues with no library coverage would silently shorten
    # the axis and shift every rolling window.
    orf_len, wt_nt, wt_aa_seq, pdb_path = 0, None, None, None
    for entry, path in zip(manifest, args.assets):
        if entry["kind"] == "aa_seq":
            wt_aa_seq = Path(path).read_text().strip().split("\\n")[0].strip()
            orf_len = len(wt_aa_seq)
        elif entry["kind"] == "fasta":
            wt_nt = wt_orf(path, run.get("reading_frame"))
        elif entry["kind"] == "pdb":
            pdb_path = path
    # A frame that disagrees with the protein length means the codon diffs would be off by a whole
    # residue, so drop the nucleotide detail rather than print something plausible and wrong.
    if wt_nt and orf_len and len(wt_nt["seq"]) // 3 not in (orf_len, orf_len + 1):
        print(f"WARNING: ORF is {len(wt_nt['seq'])} nt ({len(wt_nt['seq'])//3} codons) but the protein is "
              f"{orf_len} aa; skipping nucleotide-level detail")
        wt_nt = None

    for entry, path in zip(manifest, args.assets):
        kind, group, name = entry["kind"], entry["group"], entry["name"]

        if kind == "qc":                       # per-sample library-QC plot (original R PDF)
            # Deliberately NOT embedded. The report re-plots these natively and offers the live chart
            # as SVG/PNG, so carrying ~25 MB of base64 PDF bought nothing but weight - and the weight
            # is not free: it is what stops the embedded MultiQC report from finishing its render.
            # Only the location is kept, for the "open the results folder" links.
            plots.setdefault(group, {})[Path(name).stem] = None

        elif kind == "filtered":               # per-sample counts -> per-position curves + EC summary
            positions, classes, ec = per_position(path, orf_len)
            s = sample(group)
            s["positions"], s["classes"], s["ec"] = positions, classes, ec

        elif kind == "heatmap":                # per-sample counts-per-cov grid
            sample(group)["heatmap"] = heatmap(path)

        elif kind == "seqdepth":               # per-sample rarefaction curve (depth fraction -> % library)
            sd = pd.read_csv(path)
            sample(group)["seqdepth"] = [
                {"x": r6(x), "y": r6(y)}
                for x, y in zip(num(sd, "depth_fold"), num(sd, "variants_pct"))
                if x == x and y == y
            ]

        elif kind == "fitness_table":
            fitness = load_fitness(path, wt_nt)

        elif kind in ("fitness_pdf", "mutscan_pdf"):
            # Same reasoning as the QC PDFs: recorded by name so the page can link to them in the
            # results folder, not inlined.
            dl.append({"section": "fitness", "group": group, "name": name})

        elif kind == "report":                 # standalone detail report, embedded for a lazy iframe
            # Carried as srcdoc text rather than a data: URI. srcdoc inherits this document's origin
            # instead of minting an opaque one, and it is ~25% smaller than the base64 equivalent.
            reports.append({"group": group, "name": name, "html": html_for_frame(path)})

        elif kind == "versions":               # pipeline_info software-versions YAML
            versions = software_versions(path)

        elif kind == "trace":                  # Nextflow execution trace -> per-process compute stats
            try:
                run["stats"] = run_stats(path)
            except Exception as exc:  # a malformed/partial trace must never sink the whole report
                print(f"WARNING: could not parse execution trace: {exc}")

    # Coverage needed for every residue outcome to reach --aimed_cov counts, assuming all 21 outcomes
    # are equally represented. Same expression as the reference line in global_position_biases_cov.R;
    # it is a per-position number because each read spans the whole ORF in this assay.
    try:
        run["required_cov"] = 21 * orf_len * float(run.get("aimed_cov") or 0) or None
    except (TypeError, ValueError):
        run["required_cov"] = None
    run["orf_length"] = orf_len

    # Flags settle the tools that the version strings alone cannot: DiMSum reports only r-base, and
    # 3Dmol.js has no package entry at all - the structure viewer's presence is the tell.
    def truthy(x):
        return str(x).strip().lower() in ("true", "1", "yes")

    extra = set()
    if truthy(run.get("dimsum")):
        extra.add("dimsum")
    if truthy(run.get("mutscan")):
        extra.add("mutscan")
    if any(r["group"] == "viewer" for r in reports):
        extra.add("3dmol")
    if pdb_path:
        extra.add("3dmol")
    cites = citations(versions, extra)
    # The spinning structure needs a PDB, fitness, and the wildtype protein sequence to align them; any
    # missing and the viewer is simply left out (structure = None) so the page still renders.
    structure = structure_data(pdb_path, wt_aa_seq, fitness)
    data = {
        "run": run,
        "samples": [samples[k] for k in sorted(samples)],
        "plots": plots,
        "reports": reports,
        "fitness": fitness,
        "structure": structure,
        "downloads": dl,
        # Two citation groups for the page: the pipeline itself (always) and the tools that actually
        # ran. Each carries its own biblatex block so the "get biblatex" button has nothing to compute.
        # Section 1 is the pipeline and the frameworks it is built on; section 2 is the analysis tools
        # this particular run actually invoked.
        "citations": (lambda pipe: {
            "pipeline": [c for c in cites if c["id"] in pipe],
            "tools": [c for c in cites if c["id"] not in pipe],
            "biblatex_pipeline": biblatex([c for c in cites if c["id"] in pipe]),
            "biblatex_tools": biblatex([c for c in cites if c["id"] not in pipe]),
        })({"deepmutscan", "nfcore", "nextflow"}),
        "software": versions,
    }

    # 3Dmol.js is ~0.5 MB, so inline it only when a structure will actually be shown; otherwise the
    # placeholder collapses to nothing and a no-PDB run carries no dead weight.
    threedmol_js = ""
    if structure and args.threedmol:
        try:
            threedmol_js = Path(args.threedmol).read_text(encoding="utf-8", errors="surrogateescape")
        except OSError:
            threedmol_js = ""

    html = (Path(args.template).read_text(encoding="utf-8")
            .replace("__REPORT_DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\\\/"))
            .replace("/*__THREEDMOL_JS__*/", threedmol_js)
            .replace("__SAMPLE__", str(run.get("sample", "run"))))
    Path(args.out_html).write_text(html, encoding="utf-8")

    Path("versions.yml").write_text(
        '"${task.process}":\\n'
        + "    python: " + platform.python_version() + "\\n"
        + "    pandas: " + pd.__version__ + "\\n"
        + "    numpy: " + np.__version__ + "\\n"
    )
    mb = Path(args.out_html).stat().st_size / 1e6
    print(f"Wrote {args.out_html} ({mb:.1f} MB) - {len(samples)} samples, "
          f"{sum(len(v) for v in plots.values())} plot PDFs, {len(reports)} embedded reports")


if __name__ == "__main__":
    main()
