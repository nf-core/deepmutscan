#!/usr/bin/env python

# Nextflow template: values below are interpolated by Nextflow. The whole file is compiled
# as one Groovy GString, so every backslash must be doubled and no triple-quoted strings
# are allowed anywhere in this file, including comments.

#-- import modules --#
import os
import sys
import re
import gc
import platform
import pysam
import numpy as np
import pyarrow as pa
import polars as pl
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict, deque

#-- functions --#
def init_bam(bam_path):
    global bam_file_cov
    bam_file_cov = pysam.AlignmentFile(bam_path, "rb")

def get_base_cov(chrom, start, end, qual, min_flank):
    if min_flank and min_flank > 0:
        cov = defaultdict(int)
        for read in bam_file_cov.fetch(chrom, start, end):
            if read.is_unmapped:
                continue
            quals = read.query_qualities
            qlen = read.query_length
            if quals is None or qlen is None:
                continue
            lo = min_flank
            hi = qlen - min_flank
            for qpos, refpos in read.get_aligned_pairs(matches_only=True):
                if refpos is None or qpos is None:
                    continue
                if refpos < start or refpos >= end:
                    continue
                if qpos < lo or qpos >= hi:
                    continue
                if quals[qpos] < qual:
                    continue
                cov[refpos] += 1
        return dict(cov)

    A, C, G, T = bam_file_cov.count_coverage(chrom, start, end, quality_threshold = qual)
    coverage_array = np.array(A, dtype=np.uint32) + \\
                     np.array(C, dtype=np.uint32) + \\
                     np.array(G, dtype=np.uint32) + \\
                     np.array(T, dtype=np.uint32)

    # key on 0-based ref position to match ref_pos lookups in gatk_formating
    return {start + pos: int(cov) for pos, cov in enumerate(coverage_array)}

def chunk_ranges(start, end, sub_region_size):
    for s in range(start, end, sub_region_size):
        yield s, min(s + sub_region_size, end)

def get_base_cov_in_chunk(bam_path, chrom, start, end, qual, sub_region_size, threads, min_flank):
    chunks = list(chunk_ranges(start, end, sub_region_size))
    dicts = []

    with ProcessPoolExecutor(max_workers = threads, initializer = init_bam, initargs = (bam_path,)) as executor:
        futures = [
            executor.submit(get_base_cov, chrom, s, e, qual, min_flank)
            for s, e in chunks
        ]

        for f in as_completed(futures):
            dicts.append(f.result())

    merged_dict = {}
    for d in dicts:
        merged_dict.update(d)

    return merged_dict

# Worker state is passed explicitly through initargs rather than inherited from the parent
# address space. Relying on inheritance only works under the 'fork' start method: with 'spawn'
# (macOS) or 'forkserver' (the Python >= 3.14 default on Linux) the child re-imports this file
# with __name__ != "__main__", so names bound in the main block would be undefined.
def init_worker(base_cov, codons):
    global global_dict_base_cov, global_codon_dict
    global_dict_base_cov = base_cov
    global_codon_dict = codons

def extract_read_info(read):
    read_seq = read.query_sequence
    read_qual = read.query_qualities

    return {
        'ref':   read.reference_name,
        'pos':   read.reference_start,
        'seq':   read_seq,
        'qual':  read_qual,
        'md':    read.get_tag('MD'),
        'cigar': read.cigartuples
    }

def parse_md(read):
    read_pos = read.get('pos')
    read_md = read.get('md')
    md_pattern = re.finditer(r'(\\d+)|([A-Z]+)|(\\^[A-Z]+)', read_md)

    base_ref_pos = read_pos

    variants = []
    for match in md_pattern:
        if match.group(1):
            num_matches = int(match.group(1))
            base_ref_pos += num_matches
        elif match.group(2):
            for base in match.group(2):
                variants.append(('X', base_ref_pos, base))
                base_ref_pos += 1
        elif match.group(3):
            deleted_bases = match.group(3)[1:]
            for base in deleted_bases:
                variants.append(('D', base_ref_pos, base))
                base_ref_pos += 1

    return variants

def gatk_formating(variants, orf_start):
    if not variants:
        return 0, 0, "", 0, "", "", ""

    varying_bases = 0
    varying_codons = 0
    base_mut = ""
    codon_mut = ""
    aa_mut = ""
    pos_mut = ""

    base_pos = [ v['ref_pos'] for v in variants if v['var_type'] == 'X' ]
    base_cov = [ global_dict_base_cov.get(pos, 0) for pos in base_pos ]
    base_cov_avg = int(np.mean(base_cov)) if base_cov else 0

    base_variants = [ v for v in variants if v['var_type'] == 'X' ]
    varying_bases = len(base_variants)

    # base_mut positions are 1-based to match GATK AnalyzeSaturationMutagenesis
    base_mut = ", ".join(
        f"{v['ref_pos'] + 1}:{v['ref_base']}>{v['alt_base']}"
        for v in base_variants
    )

    codon_changes = defaultdict(list)
    for v in base_variants:
        codon_idx = ((v['ref_pos'] - orf_start) // 3) + 1
        codon_changes[codon_idx].append(v)

    varying_codons = len(codon_changes)

    codon_mut_list = []
    aa_mut_list = []
    pos_mut_list = []

    for codon_idx in sorted(codon_changes):
        vars_in_codon = codon_changes[codon_idx]
        # parse_read already resolved the reference codon for every variant it kept, keyed on the
        # read's own reference name - reuse it instead of a second lookup against a global.
        ref_codon = str(vars_in_codon[0]['ref_codon'])
        codon_list = list(ref_codon)

        for v in vars_in_codon:
            pos_in_codon = (v['ref_pos'] - orf_start) % 3
            codon_list[pos_in_codon] = v['alt_base']

        alt_codon = "".join(codon_list)
        codon_mut_list.append(f"{codon_idx}:{ref_codon}>{alt_codon}")

        ref_aa = str(Seq(ref_codon).translate())
        alt_aa = str(Seq(alt_codon).translate())

        if alt_aa == ref_aa:
            mut_type = "S"
        elif alt_aa == "*":
            mut_type = "N"
        else:
            mut_type = "M"

        # GATK writes the stop codon as X, not *
        ref_aa_s = "X" if ref_aa == "*" else ref_aa
        alt_aa_s = "X" if alt_aa == "*" else alt_aa

        aa_mut_list.append(f"{mut_type}:{ref_aa_s}>{alt_aa_s}")
        pos_mut_list.append(f"{ref_aa_s}{codon_idx}{alt_aa_s}")

    codon_mut = ", ".join(codon_mut_list)
    aa_mut = ", ".join(aa_mut_list)
    pos_mut = ";".join(pos_mut_list)

    return base_cov_avg, varying_bases, base_mut, varying_codons, codon_mut, aa_mut, pos_mut

def parse_read(read, orf_start, orf_end, base_qual, min_flank):
    read_ref = read.get('ref')
    read_pos = read.get('pos')
    read_seq = read.get('seq')
    read_qual = read.get('qual')
    read_cigar = read.get('cigar')
    read_md_variants = parse_md(read)

    if not read_md_variants:
        return []

    qlen = len(read_seq)
    base_ref_pos = read_pos
    base_read_pos = 0
    md_index = 0
    variants = []
    for op, op_len in read_cigar:
        if op == 0: # M, match
            for i in range(op_len):
                if md_index < len(read_md_variants) and read_md_variants[md_index][1] == base_ref_pos:
                    var_type, _, ref_base = read_md_variants[md_index]
                    alt_base = read_seq[base_read_pos]
                    alt_qual = read_qual[base_read_pos]
                    alt_idx = (base_ref_pos - orf_start) % 3 + 1
                    codon_idx = ((base_ref_pos - orf_start) // 3) + 1
                    variants.append({ 'var_type':  var_type,
                                      'ref_name':  read_ref,
                                      'ref_pos':   base_ref_pos,
                                      'ref_base':  ref_base,
                                      'alt_base':  alt_base,
                                      'alt_qual':  alt_qual,
                                      'alt_idx' :  alt_idx,
                                      'read_idx':  base_read_pos,
                                      'codon_idx': codon_idx })
                    md_index += 1
                base_ref_pos += 1
                base_read_pos += 1
        elif op == 1: # I, insertion (indel-containing reads are dropped upstream in bam_filter)
            base_read_pos += op_len
        elif op == 2: # D, deletion
            for i in range(op_len):
                if md_index < len(read_md_variants) and read_md_variants[md_index][1] == base_ref_pos:
                    var_type, _, ref_base = read_md_variants[md_index]
                    alt_idx = (base_ref_pos - orf_start) % 3 + 1
                    codon_idx = ((base_ref_pos - orf_start) // 3) + 1
                    variants.append({ 'var_type':  'D',
                                      'ref_name':  read_ref,
                                      'ref_pos':   base_ref_pos,
                                      'ref_base':  ref_base,
                                      'alt_base':  '-',
                                      'alt_qual':  99,
                                      'alt_idx' :  alt_idx,
                                      'read_idx':  base_read_pos,
                                      'codon_idx': codon_idx })
                    md_index += 1
                base_ref_pos += 1
        elif op == 3: # N, skip
            base_ref_pos += op_len
        elif op == 4: # S, softclip
            base_read_pos += op_len
        elif op == 5: # H, hardclip
            continue

    if not variants:
        return []

    variants_filtered = []
    for var in variants:
        if var['ref_pos'] >= orf_start and var['ref_pos'] <= orf_end:
            if var['alt_qual'] >= base_qual:
                # read-edge (min-flanking-length) filter
                if min_flank > 0:
                    read_idx = var['read_idx']
                    if read_idx < min_flank or read_idx >= qlen - min_flank:
                        continue
                codons_tmp = global_codon_dict[var['ref_name']]
                ref_codon = codons_tmp[var['codon_idx']]
                variants_filtered.append({ 'var_type':  var['var_type'],
                                           'ref_name':  var['ref_name'],
                                           'ref_pos':   var['ref_pos'],
                                           'ref_base':  var['ref_base'],
                                           'alt_base':  var['alt_base'],
                                           'alt_qual':  var['alt_qual'],
                                           'alt_idx' :  var['alt_idx'],
                                           'codon_idx': var['codon_idx'],
                                           'ref_codon': ref_codon })

    return gatk_formating(variants_filtered, orf_start)

def batch_parse_reads(batch_reads, orf_start, orf_end, base_qual, min_flank):
    cols = {
        "base_cov_avg": [],
        "varying_bases": [],
        "base_mut": [],
        "varying_codons": [],
        "codon_mut": [],
        "aa_mut": [],
        "pos_mut": []
    }

    for read in batch_reads:
        result = parse_read(read, orf_start, orf_end, base_qual, min_flank)

        if not result or result[2] == "":
            continue

        cols["base_cov_avg"].append(result[0])
        cols["varying_bases"].append(result[1])
        cols["base_mut"].append(result[2])
        cols["varying_codons"].append(result[3])
        cols["codon_mut"].append(result[4])
        cols["aa_mut"].append(result[5])
        cols["pos_mut"].append(result[6])

    if not cols["base_mut"]:
        return None

    return pa.table(cols)

def function_for_processpool(args_tuple):
    return batch_parse_reads(*args_tuple)

def aggregate_batch(table):
    # Collapse a finished worker batch to one row per distinct variant straight away. Holding one
    # row per variant-carrying read until the very end made peak memory scale with read depth;
    # aggregating here bounds it by the number of distinct variants instead.
    df = pl.from_arrow(table)
    df = df.filter(pl.col("base_mut") != "")
    if df.height == 0:
        return None
    return ( df.with_columns(pl.lit(1).alias("counts"))
               .group_by("base_mut")
               .agg([pl.col("counts").sum().alias("counts"),
                     pl.all().exclude(["base_mut", "counts"]).first()]) )

def fold_counts(acc, new):
    # Every column other than counts is a function of base_mut, so first() is value-stable here.
    if acc is None:
        return new
    return ( pl.concat([acc, new], how = "vertical")
               .group_by("base_mut")
               .agg([pl.col("counts").sum().alias("counts"),
                     pl.all().exclude(["base_mut", "counts"]).first()]) )

def read_bam_in_chunk(bam_path, orf_start, orf_end, base_qual, min_flank, batch_size, threads, base_cov, codons):
    # Keep roughly two batches per worker in flight. The previous version submitted one batch,
    # waited for it, then read the next, so only a single worker was ever busy no matter how many
    # cpus the task was given.
    max_inflight = max(2, threads * 2)

    with pysam.AlignmentFile(bam_path, "rb", threads = threads) as bam_file, \\
        ProcessPoolExecutor(max_workers = threads,
                            initializer = init_worker,
                            initargs = (base_cov, codons)) as executor:

        inflight = deque()
        batch = []

        def submit(reads):
            inflight.append(
                executor.submit(function_for_processpool,
                                (reads, orf_start, orf_end, base_qual, min_flank))
            )

        def drain_one():
            result = inflight.popleft().result()
            if result is None:
                return None
            return aggregate_batch(result)

        for read in bam_file.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.has_tag('MD'):
                continue

            # perfectly-aligned (wildtype) reads carry no variants
            if re.fullmatch(r'\\d+', read.get_tag('MD')):
                continue

            batch.append(extract_read_info(read))

            if len(batch) >= batch_size:
                submit(batch)
                batch = []
                while len(inflight) >= max_inflight:
                    agg = drain_one()
                    if agg is not None:
                        yield agg

        if batch:
            submit(batch)
            batch = []

        while inflight:
            agg = drain_one()
            if agg is not None:
                yield agg

def write_column_sheet(path):
    rows = [
        ("column", "description"),
        ("counts", "Number of reads supporting this variant (observations)"),
        ("cov", "Mean sequencing coverage across the variant's mutated base positions"),
        ("mean_length_variant_reads", "Placeholder (0); mean read length of supporting reads is not computed"),
        ("varying_bases", "Number of nucleotide substitutions in the variant"),
        ("base_mut", "Nucleotide change(s), 1-based reference position(s), e.g. 352:G>A"),
        ("varying_codons", "Number of codons altered by the variant"),
        ("codon_mut", "Codon change(s), 1-based codon index, e.g. 1:GCT>ACT"),
        ("aa_mut", "Amino-acid change per codon; type S=synonymous, M=missense, N=nonsense; stop=X, e.g. M:A>T"),
        ("pos_mut", "Compact per-codon amino-acid change, e.g. A1T (stop written as X)"),
    ]
    with open(path, "w") as fh:
        for name, desc in rows:
            fh.write(name + "\\t" + desc + "\\n")

#-- Nextflow-interpolated arguments --#
class args:
    input_bam  = "$bam"
    reference  = "$fasta"
    orf_range  = "$reading_frame"          # 1-based inclusive (GATK --orf convention)
    base_qual  = int("$base_qual")
    min_flank  = int("$min_flank")
    threads    = max(1, int("$task.cpus"))
    batch_size = 5000
    out_file   = "${meta.id}.variant_counts.tsv"

#-- main execution --#
if __name__ == "__main__":
    write_column_sheet("variant_counts_columns.tsv")

    # orf_range is 1-based inclusive (GATK --orf convention); convert to 0-based
    orf_start_1, orf_end_1 = map(int, args.orf_range.split('-'))
    orf_start = orf_start_1 - 1
    orf_end   = orf_end_1 - 1

    ref_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    orf_dict = {}
    codon_dict = {}

    if os.path.exists(args.out_file):
        os.remove(args.out_file)

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Get the reference codons, please wait...", flush = True)
    for chrom, record in ref_dict.items():
        orf_seq = record.seq[orf_start : orf_end + 1]
        orf_dict[chrom] = orf_seq
        ref_codons = {}
        for i in range(0, len(orf_seq) - 2, 3):
            codon_idx = i // 3 + 1 # 1-based codon index
            ref_codon = orf_seq[i:i + 3]
            if len(ref_codon) == 3:
                ref_codons[codon_idx] = ref_codon
        codon_dict[chrom] = ref_codons

    del ref_dict
    gc.collect()

    # The counter is built for a single gene-sized reference. Picking one sequence implicitly
    # (whatever the dict iteration happened to leave behind) would silently count against the
    # wrong contig, so require exactly one.
    if len(codon_dict) != 1:
        sys.exit(
            "ERROR: --fasta must contain exactly one sequence (found "
            + str(len(codon_dict)) + "). deepmutscan counts variants against a single ORF."
        )
    ref_name = next(iter(codon_dict))

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Get the base coverage, please wait...", flush = True)
    sub_region_size = max(1, (orf_end - orf_start + 1) // args.threads)
    dict_base_cov = get_base_cov_in_chunk(args.input_bam, ref_name, orf_start, orf_end + 1, args.base_qual, sub_region_size, args.threads, args.min_flank)

    # Nothing aligned inside the reading frame at all. Left alone this writes an empty table and
    # exits 0, and the run only dies much later in process_variant_counts.R with "no lines
    # available". Common causes: a --reading_frame outside the reference, or single-end input.
    if not dict_base_cov or sum(dict_base_cov.values()) == 0:
        sys.exit(
            "ERROR: no reads pass the quality filters inside the reading frame "
            + args.orf_range + " of " + ref_name + " for sample ${meta.id}.\\n"
            + "  Check that --reading_frame lies inside the reference, that the reads really "
            + "come from this amplicon, and that --base_qual (currently "
            + str(args.base_qual) + ") is not above the instrument quality of the data."
        )

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extracting read information, please wait...", flush = True)
    df_variants_counts = None
    for i, batch_counts in enumerate(read_bam_in_chunk(args.input_bam, orf_start, orf_end, args.base_qual, args.min_flank, args.batch_size, args.threads, dict_base_cov, codon_dict)):
        df_variants_counts = fold_counts(df_variants_counts, batch_counts)
        if (i + 1) % 100 == 0:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed batch {i+1}", flush = True)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)
    gc.collect()

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Creating the variant matrix, please wait...", flush = True)
    if df_variants_counts is None:
        df_variants_counts = pl.DataFrame([], schema={
            "base_mut": pl.Utf8,
            "counts": pl.Int64,
            "base_cov_avg": pl.Int64,
            "varying_bases": pl.Int64,
            "varying_codons": pl.Int64,
            "codon_mut": pl.Utf8,
            "aa_mut": pl.Utf8,
            "pos_mut": pl.Utf8
        }, orient = "row")

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Creating the output file, please wait...", flush = True)
    df_variants_counts = df_variants_counts.select([pl.col("counts"),
                                                    pl.col("base_cov_avg").alias("cov"),
                                                    pl.lit(0).alias("mean_length_variant_reads"),
                                                    pl.col("varying_bases"),
                                                    pl.col("base_mut"),
                                                    pl.col("varying_codons"),
                                                    pl.col("codon_mut"),
                                                    pl.col("aa_mut"),
                                                    pl.col("pos_mut")])
    # polars group_by does not promise an output order, and with several worker threads the order
    # genuinely varies between runs. That makes byte-identical reruns impossible (breaking nf-test
    # snapshots on correct data) and, worse, makes the error correction - which walks this table in
    # file order - non-reproducible. Sort on the leading reference position, then on the full
    # variant string, so the file comes out in genomic order every time.
    df_variants_counts = ( df_variants_counts
        .with_columns(pl.col("base_mut").str.extract(r"^(\\d+)", 1).cast(pl.Int64).alias("sort_pos"))
        .sort(["sort_pos", "base_mut"], nulls_last = True)
        .drop("sort_pos") )

    # --min_counts is applied at the fitness stage, not here, so the error correction
    # always sees every false double (they are typically observed only once or twice).
    # GATK variantCounts has no header line
    df_variants_counts.write_csv(args.out_file, separator = "\\t", include_header = False)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Done.", flush = True)

    #-- versions --#
    with open("versions.yml", "w") as vf:
        vf.write('"${task.process}":\\n')
        vf.write("    python: " + platform.python_version() + "\\n")
        vf.write("    pysam: " + pysam.__version__ + "\\n")
        vf.write("    polars: " + pl.__version__ + "\\n")
        vf.write("    pyarrow: " + pa.__version__ + "\\n")
        vf.write("    numpy: " + np.__version__ + "\\n")
        vf.write("    biopython: " + Bio.__version__ + "\\n")
