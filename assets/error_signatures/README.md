# Error-spectrum signature references

- `candidate_signatures_GRCh37.csv`: default candidates for `DMSANALYSIS_ERROR_SIGNATURES` (`--error_signatures_candidates`).
  One row per 96 pyrimidine-centred channel (`channel`, COSMIC notation), one column per signature, on GRCh37 trinucleotide composition.
  - SBS27, SBS43, SBS45, SBS46, SBS47, SBS50, SBS51, SBS53, SBS55, SBS57, SBS58: COSMIC v3.4 (GRCh37) signatures flagged as
    possible sequencing artefacts (Alexandrov et al. 2020, Nature 578:94; COSMIC Mutational Signatures,
    https://cancer.sanger.ac.uk/signatures/).
  - SBS18 (reactive oxygen species): COSMIC v3.4 (GRCh37). Experimentally confirmed, and oxidative damage can act on DNA in vitro.
  - FFPE_unrepaired, FFPE_repaired: formalin-fixation signatures of Guo et al. 2022 (Nat Commun 13:4487),
    `Data/FFPE_signatures_96-channel.csv` of https://github.com/QingliGuo/FFPEsig (commit 3f85d53).

  **Nine COSMIC signatures were deliberately left out of the default set**, on two grounds. Pass
  `--error_signatures_candidates assets/error_signatures/candidate_signatures_GRCh37_extended.csv` to fit with all 23.
  - _Not a process that can act on these molecules._ SBS54 is contamination with germline variants - the spectrum of common human
    polymorphisms, not a damage process. SBS30 is defective base-excision repair (NTHL1 loss), an endogenous process; it was in the
    set only as a shape proxy for unrepaired cytosine deamination, which FFPE_unrepaired represents directly (cosine 0.92).
  - _Degenerate shapes._ SBS48, SBS49, SBS52, SBS56, SBS59 and SBS60 put 40-71 % of their weight in a single trinucleotide channel
    (SBS48: 71 % in T[C>A]G; SBS60: 66 % in G[T>G]G), which is what a locus or alignment artefact looks like rather than a chemical
    process acting on every copy of a base. In a small candidate set they act as free spike absorbers.
  - SBS95 was also dropped: COSMIC reports it from a single study, with no replication across cohorts, callers or centres, and it is
    0.90 cosine-similar to SBS18, which has an experimentally confirmed mechanism. Where SBS95 was previously selected, SBS18 is now.

  Clock-like signatures (SBS1, SBS5) are deliberately not included.

- `candidate_signatures_GRCh37_extended.csv`: the full 23, including the nine above.
- `signature_descriptions.csv`: one line per signature - a short label and COSMIC's own proposed aetiology or the source publication.
  `DMSANALYSIS_ERROR_SIGNATURES` copies the entries for the signatures it selects into `error_signatures.json`, and the
  error-correction report prints them under the spectrum legend.
- `trinucleotide_counts_GRCh37.csv`: GRCh37 counts of the 32 pyrimidine-centred trinucleotides (sum over chromosomes 1-22, X, Y of
  `context_counts_GRCh37_96.csv` shipped with SigProfilerMatrixGenerator 1.4.0), used to put measured rates on the signatures' composition.

Use of COSMIC data is subject to the COSMIC terms and conditions.
