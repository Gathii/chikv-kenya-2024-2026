"""
Snakemake `script:` — summarizes the consensus-calling depth-threshold sweep.

Reads every recall_sweep/m{M}/{sample}.fa produced by the consensus_at_threshold rule,
computes real (non-N) sequence length, flags near-complete, and joins in the qPCR result 
from the provenance table. Writes both a long and a wide table.
"""
import csv

outdir = snakemake.params.outdir
thresholds = snakemake.params.thresholds
samples = snakemake.params.samples
near_complete = snakemake.params.near_complete

def real_length(fasta_path):
    seq = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    seq = "".join(seq)
    n_count = seq.upper().count("N")
    return len(seq), len(seq) - n_count


def load_pcr_status(provenance_path):
    status = {}
    with open(provenance_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sid = (row.get("Sample Id") or "").strip()
            if sid:
                status[sid.lower()] = (row.get("CHIKV qPCR Result") or "").strip()
    return status


pcr_status = load_pcr_status(snakemake.input.provenance)

long_rows = []
for sample in samples:
    lc_sample = sample.lower()
    pcr = pcr_status.get(lc_sample, "?")
    for m in thresholds:
        fa = f"{outdir}/m{m}/{sample}.fa"
        total_len, real_len = real_length(fa)
        status = "PASS" if real_len >= near_complete else "FAIL"
        long_rows.append(dict(
            sample=lc_sample, pcr_result=pcr, min_depth_threshold=m,
            total_length=total_len, real_acgt_length=real_len, status=status,
        ))

with open(snakemake.output.long, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sample", "pcr_result", "min_depth_threshold",
                                       "total_length", "real_acgt_length", "status"],
                        delimiter="\t")
    w.writeheader()
    w.writerows(long_rows)

by_sample = {}
for r in long_rows:
    by_sample.setdefault(r["sample"], {"pcr_result": r["pcr_result"]})
    mark = "" if r["status"] == "PASS" else "*"
    by_sample[r["sample"]][f"m{r['min_depth_threshold']}"] = f"{r['real_acgt_length']}{mark}"

with open(snakemake.output.wide, "w", newline="") as f:
    fieldnames = ["sample", "pcr_result"] + [f"m{m}" for m in thresholds]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    # sort by total real length across thresholds, descending (best-supported samples first)
    def sort_key(item):
        sample, row = item
        return -sum(int(row.get(f"m{m}", "0*").rstrip("*") or 0) for m in thresholds)
    for sample, row in sorted(by_sample.items(), key=sort_key):
        w.writerow({"sample": sample, **row})

print(f"Wrote {len(long_rows)} rows to {snakemake.output.long}")
print(f"Wrote {len(by_sample)} rows to {snakemake.output.wide}")
print(f"(* = below the {near_complete}bp near-complete-genome threshold)")