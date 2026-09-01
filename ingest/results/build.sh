#!/bin/bash
# Portable reproduction script for this repo's own layout (run from ingest/results/).
# Same commands/parameters as build_original_asrun.sh, which is the unmodified script as
# actually run — kept alongside for full transparency. This version just uses paths relative
# to the repo root instead of one machine's absolute paths, and expects augur/mafft/iqtree/
# treetime on PATH at the versions pinned in ../../SOFTWARE_VERSIONS.md.
set -euo pipefail
cd "$(dirname "$0")/../.."   # repo root

R=ingest/results
P=phylogenetic/results
mkdir -p "$P/alignments" "$P/tree"

echo "=== [1/9] augur align ==="
augur align \
    --method mafft \
    --sequences "$R/ecsa_seqs.fasta" \
    --reference-sequence phylogenetic/refs/CHKV.gb \
    --output "$P/alignments/ecsa_seqs-aln.fasta" \
    --fill-gaps \
    --nthreads 4

echo "=== [2/9] mask sites ==="
# sites_to_mask.txt is already provided in $P/alignments/ (same masking as the published
# build: 5' UTR 1-76, 3' UTR 11314-11826). Regenerate only if missing.
if [ ! -s "$P/alignments/sites_to_mask.txt" ]; then
    seq 1 76 > "$P/alignments/sites_to_mask.txt"
    seq 11314 11826 >> "$P/alignments/sites_to_mask.txt"
fi

augur mask \
    --sequences "$P/alignments/ecsa_seqs-aln.fasta" \
    --mask "$P/alignments/sites_to_mask.txt" \
    --output "$P/alignments/ecsa_seqs-masked.fasta"

echo "=== [3/9] augur tree (iqtree, GTR+F+I+G4, -alrt 1000) ==="
augur tree \
    --method iqtree \
    --substitution-model GTR+F+I+G4 \
    --alignment "$P/alignments/ecsa_seqs-masked.fasta" \
    --output "$P/tree/ecsa_tree_raw.nwk" \
    --tree-builder-args="-alrt 1000" \
    --nthreads 4

echo "=== [4/9] treetime clock-rate pre-check ==="
mkdir -p "$P/tree/clockcheck"
treetime clock \
    --tree "$P/tree/ecsa_tree_raw.nwk" \
    --dates "$R/metadata.tsv" \
    --aln "$P/alignments/ecsa_seqs-aln.fasta" \
    --name-column strain --date-column date \
    --outdir "$P/tree/clockcheck" 2>&1 | tee "$P/tree/clockcheck.log"

RATE=$(grep -rhoP '(?<=--rate:)\s*[0-9.eE+-]+' "$P/tree/clockcheck" "$P/tree/clockcheck.log" 2>/dev/null | tr -d ' \t' | tail -1)
if [ -z "$RATE" ]; then
    echo "WARNING: could not parse clock rate; check clockcheck.log and set manually"
    RATE=0.0003
fi
echo "Using clock rate: $RATE"
echo "$RATE" > "$P/tree/inferred_clock_rate.txt"

echo "=== [5/9] augur refine ==="
augur refine \
    --tree "$P/tree/ecsa_tree_raw.nwk" \
    --alignment "$P/alignments/ecsa_seqs-masked.fasta" \
    --metadata "$R/metadata.tsv" \
    --output-tree "$P/tree/ecsa_tree_01.nwk" \
    --output-node-data "$P/tree/branch_lengths.json" \
    --timetree \
    --coalescent opt \
    --date-confidence \
    --date-inference marginal \
    --clock-rate "$RATE" \
    --clock-std-dev 0.0001 \
    --clock-filter-iqd 2 \
    --stochastic-resolve

echo "=== [6/9] augur traits ==="
augur traits \
    --tree "$P/tree/ecsa_tree_01.nwk" \
    --metadata "$R/metadata.tsv" \
    --output-node-data "$P/traits.json" \
    --columns region country \
    --confidence \
    --sampling-bias-correction 2.5

echo "=== [7/9] augur ancestral ==="
augur ancestral \
    --tree "$P/tree/ecsa_tree_01.nwk" \
    --alignment "$P/alignments/ecsa_seqs-masked.fasta" \
    --output-node-data "$P/nt_muts.json" \
    --inference joint

echo "=== [8/9] augur translate ==="
augur translate \
    --tree "$P/tree/ecsa_tree_01.nwk" \
    --ancestral-sequences "$P/nt_muts.json" \
    --reference-sequence phylogenetic/refs/CHKV.gb \
    --output-node-data "$P/aa_muts.json" \
    --alignment-output "$P/aligned_aa_%GENE.fasta"

echo "=== [9/9] augur export v2 ==="
mkdir -p auspice
augur export v2 \
    --tree "$P/tree/ecsa_tree_01.nwk" \
    --metadata "$R/metadata.tsv" \
    --node-data "$P/tree/branch_lengths.json" "$P/traits.json" "$P/nt_muts.json" "$P/aa_muts.json" \
    --colors phylogenetic/defaults/colors.tsv \
    --lat-longs phylogenetic/defaults/lat_long.tsv \
    --auspice-config phylogenetic/defaults/auspice_config.json \
    --output auspice/ecsa_chikv.json

echo "=== BUILD COMPLETE: auspice/ecsa_chikv.json ==="
