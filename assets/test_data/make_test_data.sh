#!/usr/bin/env bash
# Generates throwaway inputs for the -stub-run CI tests.
# Nothing here is real data: the FASTQs are empty gzip streams and the
# "reference" is an empty directory. Only the channel logic is exercised.

set -euo pipefail

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FQ="${HERE}/fastq"
REF="${HERE}/ref"

mkdir -p "$FQ" "$REF"

# fake transcriptome dir - cellranger never runs in stub mode
mkdir -p "$REF/star" "$REF/fasta" "$REF/genes"
echo '{"genomes": ["testgenome"], "mem_gb": 1}' > "$REF/reference.json"

make_fq () {
    printf '' | gzip -c > "${FQ}/$1"
}

# sample A: one library, one run
make_fq sampleA_run1_R1.fastq.gz
make_fq sampleA_run1_R2.fastq.gz

# sample B: one library, TWO sequencing runs (top-up) -> merged as pseudo-lanes
make_fq sampleB_run1_R1.fastq.gz
make_fq sampleB_run1_R2.fastq.gz
make_fq sampleB_run2_R1.fastq.gz
make_fq sampleB_run2_R2.fastq.gz

# sample C: TWO libraries (separate GEM wells) -> separate counts, same sample
make_fq sampleC_libA_R1.fastq.gz
make_fq sampleC_libA_R2.fastq.gz
make_fq sampleC_libB_R1.fastq.gz
make_fq sampleC_libB_R2.fastq.gz

# multiome sample: gex pair + atac triple
make_fq brain1_gex_R1.fastq.gz
make_fq brain1_gex_R2.fastq.gz
make_fq brain1_atac_R1.fastq.gz
make_fq brain1_atac_R2.fastq.gz
make_fq brain1_atac_R3.fastq.gz

# ---- samplesheets, written with absolute paths -----------------------------

cat > "${HERE}/samplesheet_gex.csv" <<EOF
sample,library,fastq_1,fastq_2,modality
sampleA,,${FQ}/sampleA_run1_R1.fastq.gz,${FQ}/sampleA_run1_R2.fastq.gz,gex
EOF

cat > "${HERE}/samplesheet_multi.csv" <<EOF
sample,library,fastq_1,fastq_2,modality
sampleA,,${FQ}/sampleA_run1_R1.fastq.gz,${FQ}/sampleA_run1_R2.fastq.gz,gex
sampleB,,${FQ}/sampleB_run1_R1.fastq.gz,${FQ}/sampleB_run1_R2.fastq.gz,gex
sampleB,,${FQ}/sampleB_run2_R1.fastq.gz,${FQ}/sampleB_run2_R2.fastq.gz,gex
sampleC,sampleC_libA,${FQ}/sampleC_libA_R1.fastq.gz,${FQ}/sampleC_libA_R2.fastq.gz,gex
sampleC,sampleC_libB,${FQ}/sampleC_libB_R1.fastq.gz,${FQ}/sampleC_libB_R2.fastq.gz,gex
EOF

cat > "${HERE}/samplesheet_arc.csv" <<EOF
sample,library,fastq_1,fastq_2,fastq_3,modality
brain1,brain1_gex,${FQ}/brain1_gex_R1.fastq.gz,${FQ}/brain1_gex_R2.fastq.gz,,arc_gex
brain1,brain1_atac,${FQ}/brain1_atac_R1.fastq.gz,${FQ}/brain1_atac_R2.fastq.gz,${FQ}/brain1_atac_R3.fastq.gz,arc_atac
EOF

echo "Test data written to ${HERE}"
