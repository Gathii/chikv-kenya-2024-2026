configfile: "config.yaml"

# Define final outputs
rule all:
    input:
        f"auspice/{config['build']['name']}.json"

# 1. Indexing & Filtering
# -----------------------------------------------------------------------------

rule index_context:
    input:
        sequences = config["inputs"]["context_seqs"]
    output:
        index = "ingest/results/context_seqs_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.index}
        """

rule filter_south_america:
    input:
        sequences = config["inputs"]["context_seqs"],
        metadata = config["inputs"]["context_meta"],
        index = "ingest/results/context_seqs_index.tsv"
    output:
        sequences = "ingest/results/subsampled_south_america.fasta"
    params:
        filter_str = config["filtering"]["south_america"]["query"],
        group_by = config["filtering"]["south_america"]["group_by"],
        per_group = config["filtering"]["south_america"]["seqs_per_group"],
        min_len = config["filtering"]["min_length"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude-where {params.filter_str} \
            --group-by {params.group_by} \
            --sequences-per-group {params.per_group} \
            --min-length {params.min_len} \
            --output-sequences {output.sequences}
        """

rule filter_kenya:
    input:
        sequences = config["inputs"]["context_seqs"],
        metadata = config["inputs"]["context_meta"],
        index = "ingest/results/context_seqs_index.tsv"
    output:
        sequences = "ingest/results/subsampled_kenya.fasta"
    params:
        filter_str = config["filtering"]["kenya"]["query"],
        group_by = config["filtering"]["kenya"]["group_by"],
        per_group = config["filtering"]["kenya"]["seqs_per_group"],
        min_len = config["filtering"]["kenya"]["min_length"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude-where {params.filter_str} \
            --group-by {params.group_by} \
            --sequences-per-group {params.per_group} \
            --min-length {params.min_len} \
            --output-sequences {output.sequences}
        """

rule filter_global:
    input:
        sequences = config["inputs"]["context_seqs"],
        metadata = config["inputs"]["context_meta"],
        index = "ingest/results/context_seqs_index.tsv"
    output:
        sequences = "ingest/results/global_background.fasta"
    params:
        filter_str = config["filtering"]["global"]["query"],
        group_by = config["filtering"]["global"]["group_by"],
        per_group = config["filtering"]["global"]["seqs_per_group"],
        min_len = config["filtering"]["min_length"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude-where {params.filter_str} \
            --group-by {params.group_by} \
            --sequences-per-group {params.per_group} \
            --min-length {params.min_len} \
            --output-sequences {output.sequences}
        """

rule index_focal:
    input:
        sequences = config["inputs"]["focal_seqs"]
    output:
        index = "ingest/results/focal_seqs_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.index}
        """

rule filter_focal:
    input:
        sequences = config["inputs"]["focal_seqs"],
        metadata = config["inputs"]["focal_meta"],
        index = "ingest/results/focal_seqs_index.tsv"
    output:
        sequences = "ingest/results/focal_filt_seqs.fasta"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --output-sequences {output.sequences} \
            --min-date 2000 \
            --min-length 8000
        """

# 2. Combine Data
# -----------------------------------------------------------------------------

rule combine_sequences:
    input:
        sa = "ingest/results/subsampled_south_america.fasta",
        ke = "ingest/results/subsampled_kenya.fasta",
        gl = "ingest/results/global_background.fasta",
        fo = "ingest/results/focal_filt_seqs.fasta"
    output:
        "ingest/results/ecsa_seqs.fasta"
    shell:
        "cat {input.sa} {input.ke} {input.gl} {input.fo} > {output}"

rule combine_metadata:
    input:
        focal = config["inputs"]["focal_meta"],
        context = config["inputs"]["context_meta"]
    output:
        "ingest/results/metadata.tsv"
    shell:
        # Warning: Simple cat requires identical headers in both files
        "cat {input.focal} {input.context} > {output}"

# 3. Alignment & Masking
# -----------------------------------------------------------------------------

rule align:
    input:
        sequences = "ingest/results/ecsa_seqs.fasta",
        ref = config["inputs"]["reference"]
    output:
        alignment = "phylogenetic/results/alignments/ecsa_seqs-aln.fasta"
    threads: 16
    shell:
        """
        augur align \
            --method mafft \
            --sequences {input.sequences} \
            --reference-sequence {input.ref} \
            --output {output.alignment} \
            --fill-gaps \
            --nthreads {threads}
        """

rule generate_mask_file:
    output:
        mask = "phylogenetic/results/alignments/sites_to_mask.txt"
    shell:
        """
        seq 1 76 > {output.mask} && \
        seq 11314 11826 >> {output.mask}
        """

rule mask:
    input:
        alignment = "phylogenetic/results/alignments/ecsa_seqs-aln.fasta",
        mask_sites = "phylogenetic/results/alignments/sites_to_mask.txt"
    output:
        masked = "phylogenetic/results/alignments/ecsa_seqs-masked.fasta"
    shell:
        """
        augur mask \
            --sequences {input.alignment} \
            --mask {input.mask_sites} \
            --output {output.masked}
        """

# 4. Tree Building
# -----------------------------------------------------------------------------

rule tree:
    input:
        alignment = "phylogenetic/results/alignments/ecsa_seqs-masked.fasta"
    output:
        tree = "phylogenetic/results/tree/ecsa_tree_raw.nwk"
    threads: 20
    shell:
        """
        augur tree \
            --method iqtree \
            --substitution-model GTR+F+I+G4 \
            --alignment {input.alignment} \
            --output {output.tree} \
            --tree-builder-args="-alrt 1000" \
            --nthreads {threads}
        """

rule clock_check:
    input:
        tree = "phylogenetic/results/tree/ecsa_tree_raw.nwk",
        alignment = "phylogenetic/results/alignments/ecsa_seqs-aln.fasta",
        metadata = "ingest/results/metadata.tsv"
    output:
        rate_file = "phylogenetic/results/tree/inferred_clock_rate.txt"
    params:
        outdir = "phylogenetic/results/tree/clockcheck",
        fallback = config["refine"]["clock_rate_fallback"]
    shell:
        """
        mkdir -p {params.outdir}
        treetime clock \
            --tree {input.tree} \
            --dates {input.metadata} \
            --aln {input.alignment} \
            --name-column strain --date-column date \
            --outdir {params.outdir} 2>&1 | tee phylogenetic/results/tree/clockcheck.log
        RATE=$(grep -rhoP '(?<=--rate:)\\s*[0-9.eE+-]+' {params.outdir} phylogenetic/results/tree/clockcheck.log 2>/dev/null | tr -d ' \\t' | tail -1)
        if [ -z "$RATE" ]; then
            echo "WARNING: could not parse clock rate from treetime output; using config fallback {params.fallback}" >&2
            RATE={params.fallback}
        fi
        echo "$RATE" > {output.rate_file}
        """

rule refine:
    input:
        tree = "phylogenetic/results/tree/ecsa_tree_raw.nwk",
        alignment = "phylogenetic/results/alignments/ecsa_seqs-masked.fasta",
        metadata = "ingest/results/metadata.tsv",
        rate_file = "phylogenetic/results/tree/inferred_clock_rate.txt"
    output:
        tree = "phylogenetic/results/tree/ecsa_tree_01.nwk",
        node_data = "phylogenetic/results/tree/branch_lengths.json"
    params:
        std_dev = config["refine"]["clock_std_dev"]
    shell:
        """
        RATE=$(cat {input.rate_file})
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent opt \
            --date-confidence \
            --date-inference marginal \
            --clock-rate "$RATE" \
            --clock-std-dev {params.std_dev} \
            --clock-filter-iqd 2 \
            --stochastic-resolve
        """

# 5. Annotation & Export
# -----------------------------------------------------------------------------

rule traits:
    input:
        tree = "phylogenetic/results/tree/ecsa_tree_01.nwk",
        metadata = "ingest/results/metadata.tsv"
    output:
        node_data = "phylogenetic/results/traits.json"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns region country \
            --confidence \
            --sampling-bias-correction 2.5
        """

rule ancestral:
    input:
        tree = "phylogenetic/results/tree/ecsa_tree_01.nwk",
        alignment = "phylogenetic/results/alignments/ecsa_seqs-masked.fasta"
    output:
        node_data = "phylogenetic/results/nt_muts.json"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference joint
        """

rule translate:
    input:
        tree = "phylogenetic/results/tree/ecsa_tree_01.nwk",
        ancestral = "phylogenetic/results/nt_muts.json",
        ref = config["inputs"]["reference"]
    output:
        node_data = "phylogenetic/results/aa_muts.json",
        # Using a marker file because 'aligned_aa_%GENE.fasta' yields multiple files
        marker = "phylogenetic/results/.aa_translation_done"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.ancestral} \
            --reference-sequence {input.ref} \
            --output-node-data {output.node_data} \
            --alignment-output phylogenetic/results/aligned_aa_%GENE.fasta && \
        touch {output.marker}
        """

rule export:
    input:
        tree = "phylogenetic/results/tree/ecsa_tree_01.nwk",
        metadata = "ingest/results/metadata.tsv",
        branch_lengths = "phylogenetic/results/tree/branch_lengths.json",
        traits = "phylogenetic/results/traits.json",
        nt_muts = "phylogenetic/results/nt_muts.json",
        aa_muts = "phylogenetic/results/aa_muts.json",
        colors = config["inputs"]["colors"],
        lat_longs = config["inputs"]["lat_longs"],
        auspice_config = config["inputs"]["auspice_config"]
    output:
        auspice_json = f"auspice/{config['build']['name']}.json"
    params:
        title = config["build"]["title"],
        maintainers = config["build"]["maintainers"],
        url = config["build"]["url"]
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} \
                        {input.traits} \
                        {input.nt_muts} \
                        {input.aa_muts} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --metadata-columns country region division location year lineage \
            --output {output.auspice_json} \
            --title "{params.title}" \
            --maintainers "{params.maintainers}" \
            --build-url "{params.url}"
        """
