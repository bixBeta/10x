#!/usr/bin/env bash
# Throwaway inputs for the -stub-run CI check. The fastqs are empty gzip
# streams and the reference is an empty dir; only the channel logic is tested.

set -euo pipefail

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FQ="${HERE}/fastq"
REF="${HERE}/ref"

mkdir -p "$FQ" "$REF"
echo '{"genomes": ["testgenome"]}' > "$REF/reference.json"

mk () { printf '' | gzip -c > "${FQ}/$1" ; }

# SS1: one library, one run
mk SS1_R1.fastq.gz ; mk SS1_R2.fastq.gz
# SS2: one library, two sequencing runs -> pooled as L001 + L002
mk SS2_run1_R1.fastq.gz ; mk SS2_run1_R2.fastq.gz
mk SS2_run2_R1.fastq.gz ; mk SS2_run2_R2.fastq.gz
# SS3: two libraries / GEM wells -> two separate counts
mk SS3_A_R1.fastq.gz ; mk SS3_A_R2.fastq.gz
mk SS3_B_R1.fastq.gz ; mk SS3_B_R2.fastq.gz

cat > "${HERE}/sample-sheet.csv" <<EOF
label,library,fastq1,fastq2
SS1,,${FQ}/SS1_R1.fastq.gz,${FQ}/SS1_R2.fastq.gz
SS2,,${FQ}/SS2_run1_R1.fastq.gz,${FQ}/SS2_run1_R2.fastq.gz
SS2,,${FQ}/SS2_run2_R1.fastq.gz,${FQ}/SS2_run2_R2.fastq.gz
SS3,SS3_wellA,${FQ}/SS3_A_R1.fastq.gz,${FQ}/SS3_A_R2.fastq.gz
SS3,SS3_wellB,${FQ}/SS3_B_R1.fastq.gz,${FQ}/SS3_B_R2.fastq.gz
EOF

echo "test data written to ${HERE}"
