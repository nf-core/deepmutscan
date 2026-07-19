# nf-core/deepmutscan — all-in-one report: wishes & status

The report is one self-contained HTML file (`deepmutscan_report.html`), built by
`VISUALIZATION_SUMMARY_REPORT`, covering QC → fitness for every sequencing file; intuitive,
interactive, modern/minimalist, scientifically correct, not "AI-generated"-looking. Heading styled
like the inspection tool. Info buttons throughout.

## Done & verified

- Self-contained single file; embeds MultiQC + EC reports; re-plots tables interactively.
- MultiQC embed renders (srcdoc + localStorage/cookie shims for the opaque-origin frame).
- Light mode default, dark optional via header toggle; dual logo (white-font logo in dark).
- Theme propagates into embedded MultiQC + EC reports.
- Axis titles no longer overlap tick labels (margins measured, not hardcoded).
- Library QC curves split by varying-bases class (1/2/3) matching the pipeline's R semantics
  (within-class means, rolling mean over `--sliding_window_size` on by default + toggle,
  required-coverage reference line). Class selector is multi-select with transparency: default
  "1 varying base" opaque, others faded; any combination selectable.
- Per-plot downloads: high-res PNG + self-contained (dependency-free) PDF; no SVG. Word "original"
  removed everywhere (report plots are original too).
- Downloads page removed; results-folder links removed; Overview has a README (purpose + output
  file-structure table).
- "Mean fitness along the ORF" removed.
- Fitness page renamed **Variant effects**.
- Fitness landscape heatmap is interactive: click a cell → exact nucleotide change(s) contributing
  (validated 10,254/10,254 against the pipeline's `base_mut`) + per-replicate fitness, SD, in/out
  counts.
- Citation page: section 1 pipeline (Wehnert et al., bioRxiv preprint coming soon) + nf-core +
  Nextflow; section 2 tools actually used this run; each with a "Get biblatex" button opening a new
  tab. (Also fixed 3 modules emitting a literal `${task.process}` into `versions.yml`, and corrected
  the BWA citation to BWA-MEM arXiv:1303.3997.)
- **EC report redesign — DONE** (no longer pools). Three-section structure:
  1. **Compare** — a per-file correlation matrix of the sequencing-error signature: for each SNV the
     error frequency removed (raw − corrected counts-per-coverage, which equals the per-nt error rate
     since coverage cancels — so no extra file has to be threaded through the pipeline) is one number
     per file, and the matrix is the Pearson r of those between every file pair. High r = the same
     positions/substitutions are error-prone in every file (systematic chemistry). Works for both
     false_doubles and wildtype; degrades to a "needs ≥2 files" note.
  2. **Selection** — "Raw vs corrected frequency" split into an **input-libraries** plot and an
     **output-libraries** plot, each with its own file toggles; grouping is by the samplesheet
     `type`/`replicate` (general, not GID1A-specific — threaded via the report's file metadata).
     Multi-file selections draw a spread bar.
  3. **Detail** — a single-file segmented selector drives: positional error bias, a per-class
     beeswarm of raw (filled) + corrected (hollow) frequencies ("Sequencing-error frequency by
     substitution class"), a **correction-magnitude** panel now plotting the *frequency removed* per
     adjusted variant on a log axis with fine bins (the old %-change was bimodal at 0 / −100 % because
     the estimator either leaves a variant alone or floors it to zero), and the SNV table.
  Tooltips are delegated (one handler per SVG) to keep ~4 k points light. Verified against
  results_GID1A via the harness; workflow compiles (`-preview`, exit 0). One wiring change only:
  the report now carries `{id,type,replicate}` per file instead of just the id.

## Done since (this session)

- **Report opens on the Variant effects page** when the run has fitness (falls back to Overview
  otherwise) — DONE.
- **Variant effects: results-folder note** — DONE. A "Complete fitness outputs" panel lists each
  estimator that ran and where its full tables/plots/report live. (Spinning 3Dmol canvas + inspection
  link still outstanding, below.)
- **EC report refinements** — DONE: correction-magnitude x-axis is now "% removed" (0–100, linear);
  the "by substitution class" beeswarm was reordered above the ORF-bias panel; section-3 heading + file
  selector are now sticky while its panels scroll; the beeswarm no longer overlaps (single full-width
  swarm per class + tie-breaking jitter, no compression — 2 grazing pairs vs 628 before); all info
  buttons trimmed to ≤2 sentences describing only what each plot is built from.
  - NOTE on "% removed": the false-doubles estimator is all-or-nothing on this data (a variant is
    either untouched or removed entirely), so the histogram is essentially one bar at 100 %. Honest but
    visually flat; the alternative "frequency removed (log)" view is more graded if a change of mind.

- **Variant effects: spinning protein** — DONE. A live 3Dmol canvas (transparent background, gentle
  auto-spin) shows the wildtype fold coloured by mean fitness per residue, aligned to the sequence
  (identity shown when <90%); plus an "Open the full interactive inspection tool ↗" link that opens the
  embedded viewer from a blob, so it survives the report being downloaded/sent. 3Dmol.js is inlined by
  the build only when the run has a `--pdb` (no dead weight otherwise). Wiring: PDB threaded into the
  report manifest; `3Dmol-min.js` passed to VISUALIZATION_SUMMARY_REPORT.
- **Sequencing QC: seqdepth own-data plot** — DONE. `SeqDepth_simulation.R` rewritten to closed-form
  hypergeometric rarefaction (verified byte-identical to an independent numpy implementation), emits a
  curve CSV, `run_seqdepth` now defaults **true**; the report re-plots all files as overlaid saturation
  curves under Sequencing QC.
- **Overview run statistics** — DONE (no CO2, per Benjamin). Parses the Nextflow execution trace into
  per-process tasks / cached / CPU time / mean %CPU / peak RSS + summary tiles; labelled "as of report
  generation".

## Done (round 2, 2026-07-18)

- **Report published at the results root** (`deepmutscan_report.html`), no longer under
  `library_QC/run/`. Needed both a new `withName: 'VISUALIZATION_SUMMARY_REPORT'` publishDir *and*
  excluding it from the generic `VISUALIZATION_*` selector — Nextflow accumulates publishDir across
  matching selectors rather than overriding.
- **Opens on Overview** again (was Variant effects).
- **Spinning protein**: panel 440 → 820 px; canvas `pointer-events: none` (interaction lives in the
  inspection tool); and the viewer is now **built once at report load into an off-screen holder** and
  moved into the page on show / rescued back before each page switch — so it is already warmed up and
  spinning on arrival and its WebGL context is never rebuilt.
- **EC report loads faster**: the heavy single-file detail section (bias, beeswarm, magnitude, table)
  is deferred to the next tick, so the comparison panels paint first.
- **Beeswarm side-by-side**: raw left / corrected right, both filled. Placement is now a **single pass
  over every class and lane** (absolute x) — placing lanes or classes separately let dense swarms bleed
  into their neighbours. Verified **0 overlapping points** (min separation 2.7 px vs 2.2 px markers).
- **Correction magnitude** now includes untouched variants, so the 0 % bar is there (917 on input_1).
- **Library QC** heading + sequencing-file selector are sticky while the plots scroll.
- **Logo** 34 → 50 px. **Downloads are PNG only** (PDF button and the PDF builder removed).

## Done (round 3, 2026-07-19)

- **Correction magnitude shows corrected variants only** again (the 0 % bar is gone), and its info button
  now says so explicitly. NOTE: with untouched variants excluded the panel is a single bar at 98-100 %
  on this data, because the false-doubles estimator is all-or-nothing (a variant is either left alone or
  floored to zero).
- **Every page opens at the top.** `main` is the scroll container, so its `scrollTop` survived a page
  switch and pages opened at whatever offset the previous one had been left at; `show()` now resets it.
- **Spinning protein is smooth.** The real cause was not the GPU: 3Dmol's `viewer.spin()` drives rotation
  from `setInterval(..., 25)` - a hard 40 fps ceiling that is also not vsync-aligned, so on a 60 Hz display
  frames land one or two refreshes apart and visibly judder. It is now driven by `requestAnimationFrame`
  (vsync-locked, time-based so the speed is identical on any refresh rate, step clamped so a backgrounded
  tab does not snap-rotate on return) and paused while the panel is off-screen. Also `antialias: false` +
  `cartoonQuality: 4` (from 10): 3Dmol's antialias additionally *doubles* the internal devicePixelRatio.
- **Whitespace above/below the structure** removed: panel 820 -> 560 px and `zoom(1.35)` after `zoomTo()`,
  which fits the whole model and so left the fold floating in dead space on a wide panel.

## Outstanding

- **Method resilience:** largely by construction (every page/section is guarded — mqc, fitness,
  structure, EC, each estimator); a from-scratch no-fitness / no-PDB / single-method build still worth
  one explicit verification pass.
- **Titles/subtitles pass** across all pages — light touch done (Overview lead strengthened, EC titles
  reworked, new panels titled); a full editorial sweep is still open.
- **rollmean vs R** still validated only by reading rollapply's docs, not by testing — verify against R
  before fully trusting the smoothed Library QC curves.
- **CO2** (nf-co2footprint) intentionally deferred by Benjamin (2026-07-18).
- **Secondary-structure track** — DROPPED per Benjamin (2026-07-18).

## Containers needed for outstanding report features

One shared image: `pandas=2.2.1`, `numpy` (let float; biotite may need ≥2.0), `biotite` (secondary
structure + inspection tool), `pymupdf` (thumbnail rasterising).
