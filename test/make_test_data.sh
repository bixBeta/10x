#!/usr/bin/env bash
# Throwaway inputs for CI. Builds 10x delivery dirs with the real naming
# convention so the --sample derivation is exercised for what it is:
#
#   Sample_<prefix>_L<lane>/<prefix>_S<n>_L00<lane>_R{1,2}_001.fastq.gz
#
# The fastqs are empty gzip streams; only the naming matters.

set -euo pipefail

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RUN1="${HERE}/Illumina/run1/Unaligned/Project_10488923"
RUN2="${HERE}/Illumina/run2/Unaligned/Project_10488923"
REF="${HERE}/ref"

rm -rf "${HERE}/Illumina"
mkdir -p "$RUN1" "$RUN2" "$REF"
echo '{"genomes": ["testgenome"]}' > "$REF/reference.json"

# make_delivery <parent> <prefix> <S-index> <lane>
make_delivery () {
    local parent="$1" prefix="$2" sidx="$3" lane="$4"
    local d="${parent}/Sample_${prefix}_L${lane}"
    mkdir -p "$d"
    for r in R1 R2 ; do
        printf '' | gzip -c > "${d}/${prefix}_S${sidx}_L00${lane}_${r}_001.fastq.gz"
        : > "${d}/${prefix}_S${sidx}_L00${lane}_${r}_001.md5"
    done
    : > "${d}/${prefix}_S${sidx}_L00${lane}.fastp.json"
    echo "$d"
}

# JS4: one library, one delivery dir
JS4=$( make_delivery "$RUN1" "SC2620_JS4_G3_Reign_10488923_253GGLLT4" 3 2 )

# JS5: one library, two sequencing runs. The flow cell is part of the prefix,
# so the two runs have DIFFERENT sample prefixes and both must be listed.
JS5a=$( make_delivery "$RUN1" "SC2620_JS5_G3_Sofia_10488923_253GGLLT4" 4 2 )
JS5b=$( make_delivery "$RUN2" "SC2620_JS5_G3_Sofia_10488923_999XYZAB2" 1 3 )

# JS6: two libraries / GEM wells under one label
JS6a=$( make_delivery "$RUN1" "SC2620_JS6A_G3_Scooter_10488923_253GGLLT4" 5 2 )
JS6b=$( make_delivery "$RUN1" "SC2620_JS6B_G3_Scooter_10488923_253GGLLT4" 6 2 )

# multiome: one label, a GEX library and an ATAC library, delivered under one
# Project dir - which is what the real libraries.csv points --fastqs at
ARCP="${HERE}/Illumina/arcrun/Unaligned/Project_10488522"
mkdir -p "$ARCP"
make_delivery "$ARCP" "SC2619_JS4_BC_MG3_8Healthy_10488522_25FWVCLT4" 1 1 > /dev/null
make_delivery "$ARCP" "SC2619_JS4_MA_8Healthy_10488522_23C52HLT4"     2 1 > /dev/null

# the Project dir holds both libraries, so each row names its sample - exactly
# as the hand written libraries.csv does
cat > "${HERE}/sample-sheet-arc.csv" <<EOF
label,fastqs,sample,library_type
JS4,${ARCP},SC2619_JS4_BC_MG3_8Healthy_10488522_25FWVCLT4,Gene Expression
JS4,${ARCP},SC2619_JS4_MA_8Healthy_10488522_23C52HLT4,Chromatin Accessibility
EOF

# same Project dir with no sample column: ambiguous, must be rejected
cat > "${HERE}/sample-sheet-ambiguous.csv" <<EOF
label,fastqs
JS4,${ARCP}
EOF

cat > "${HERE}/sample-sheet.csv" <<EOF
label,library,fastqs
JS4,,${JS4}
JS5,,${JS5a}
JS5,,${JS5b}
JS6,JS6_wellA,${JS6a}
JS6,JS6_wellB,${JS6b}
EOF

echo "test data written to ${HERE}"
