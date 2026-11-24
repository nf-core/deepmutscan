#!/usr/bin/env bash
#
# tile_flanks.sh REFERENCE_LENGTH < SAM records
#
# nf-core/deepmutscan, --sequencing tiled: finds the constant sequence that the merged reads of a tiled
# amplicon carry beyond each end of the reference, from a sample of their local alignments (bwa-mem3
# default scoring, which soft-clips what the reference does not contain). Input: SAM records without
# header, primary alignments only.
#
# At the 5' end, the clipped sequences of the reads whose alignment starts at reference base 1 are read
# backwards from the reference boundary, at the 3' end those of the reads whose alignment ends at the
# last reference base forwards from it. The flank grows base by base for as long as at least half of
# those reads are still long enough and >= 90 % of them carry the same base; a variable spacer or
# barcode beyond the constant flank stops it. An end is only considered when at least half of all
# aligned reads are clipped there.
#
# Output, tab-separated, one line per end: end (5prime|3prime), flank (empty if none), reads clipped
# at that end, reads aligned.

set -euo pipefail

awk -v len="$1" -v min_share=0.5 -v min_major=0.9 '
    # one base of the flank at distance j from the reference boundary, or "" to stop
    function base_at(clips, n, j, five,    i, c, b, tot, cnt, best, bmax) {
        tot = 0; split("", cnt)
        for (i = 1; i <= n; i++) {
            c = clips[i]
            if (length(c) < j) continue
            b = five ? substr(c, length(c) - j + 1, 1) : substr(c, j, 1)
            cnt[b]++; tot++
        }
        if (tot == 0 || tot < min_share * n) return ""
        bmax = 0; best = ""
        for (b in cnt) if (cnt[b] > bmax) { bmax = cnt[b]; best = b }
        return (bmax >= min_major * tot) ? best : ""
    }
    function flank(clips, n, five,    j, b, f) {
        f = ""
        if (n < min_share * aligned) return f
        for (j = 1; ; j++) {
            b = base_at(clips, n, j, five)
            if (b == "" || b == "N") break
            f = five ? b f : f b
        }
        return f
    }
    {
        aligned++
        cigar = $6; seq = $10
        if ($4 == 1 && match(cigar, /^[0-9]+S/)) c5[++n5] = substr(seq, 1, substr(cigar, 1, RLENGTH - 1) + 0)
        if (match(cigar, /[0-9]+S$/)) {
            clip = substr(cigar, RSTART, RLENGTH - 1) + 0
            end = $4 - 1; t = cigar
            while (match(t, /^[0-9]+[MIDNSHP=X]/)) {
                op = substr(t, RLENGTH, 1)
                if (op ~ /[MDN=X]/) end += substr(t, 1, RLENGTH - 1) + 0
                t = substr(t, RLENGTH + 1)
            }
            if (end == len) c3[++n3] = substr(seq, length(seq) - clip + 1)
        }
    }
    END {
        printf "5prime\t%s\t%d\t%d\n", flank(c5, n5, 1), n5, aligned
        printf "3prime\t%s\t%d\t%d\n", flank(c3, n3, 0), n3, aligned
    }
'
