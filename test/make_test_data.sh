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
rm -rf "${HERE}/refparent"

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

# scATAC: one row per library, same shape as gex. The second label is a
# re-sequenced library, so --fastqs / --sample have to come out as lists.
cat > "${HERE}/sample-sheet-atac.csv" <<EOF
label,fastqs,library_type
CK1,${JS4},Chromatin Accessibility
CK2,${JS5a},Chromatin Accessibility
CK2,${JS5b},Chromatin Accessibility
EOF

# multi: one label carrying GEX + Antibody Capture + VDJ, the shape of a real
# cellranger multi submission. All three live in one Project dir, so each row
# names its fastq_id.
MULTIP="${HERE}/Illumina/multirun/Unaligned/Project_10489440"
mkdir -p "$MULTIP"
make_delivery "$MULTIP" "SC2622_EW1_G5_HL978_10489440_25GHMFLT4"  1 2 > /dev/null
make_delivery "$MULTIP" "SC2622_EW1_CSP_HL978_10489440_25GHMFLT4" 2 2 > /dev/null
make_delivery "$MULTIP" "SC2622_EW1_VDJ_HL978_10489440_25GHMFLT4" 3 5 > /dev/null

# the feature barcode panel the [feature] section points at
cat > "${HERE}/feature_ref.csv" <<EOF
id,name,read,pattern,sequence,feature_type
CD3,CD3,R2,^NNNNNNNNNN(BC)NNNNNNNNN,AACAAGACCCTTGAG,Antibody Capture
EOF

# a vdj reference, checked the same way a gex one is
mkdir -p "${HERE}/vdjref"
echo '{"genomes": ["testvdj"]}' > "${HERE}/vdjref/reference.json"

cat > "${HERE}/sample-sheet-multi.csv" <<EOF
label,fastqs,sample,lanes,library_type
EW1,${MULTIP},SC2622_EW1_G5_HL978_10489440_25GHMFLT4,2,Gene Expression
EW1,${MULTIP},SC2622_EW1_CSP_HL978_10489440_25GHMFLT4,2,Antibody Capture
EW1,${MULTIP},SC2622_EW1_VDJ_HL978_10489440_25GHMFLT4,5,VDJ
EOF

# the same sheet without a VDJ reference available: must be rejected
cat > "${HERE}/sample-sheet-multi-novdjref.csv" <<EOF
label,fastqs,sample,library_type
EW1,${MULTIP},SC2622_EW1_G5_HL978_10489440_25GHMFLT4,Gene Expression
EW1,${MULTIP},SC2622_EW1_VDJ_HL978_10489440_25GHMFLT4,VDJ
EOF

# the sample column filled in with a file name stem, illumina tail and all -
# the error should name the prefix meant rather than only report nothing found
cat > "${HERE}/sample-sheet-filestem.csv" <<EOF
label,fastqs,sample
JS4,${JS4},SC2620_JS4_G3_Reign_10488923_253GGLLT4_S3_L002
EOF

# visium: one capture area per label, delivered one directory PER LANE, which is
# what the real spatial deliveries look like. YF1 spans two lanes, so its rows
# repeat and have to pool; YF2 is a second capture area on the same slide.
SPATP="${HERE}/Illumina/spatialrun/Unaligned/Project_10490358"
mkdir -p "$SPATP"
YF1a=$( make_delivery "$SPATP" "SC2630_YF1_V-G3_A1_10490358_23NY35LT3" 10 1 )
YF1b=$( make_delivery "$SPATP" "SC2630_YF1_V-G3_A1_10490358_23NY35LT3" 10 2 )
YF2=$(  make_delivery "$SPATP" "SC2630_YF2_V-G3_B1_10490358_23NY35LT3" 11 1 )

# the images the lab delivers alongside the reads, dropped in one folder. The
# sheet names them, --imagedir says where to look.
# CytAssist writes <instrument>_<date>_<time>_<date>_<time>_<slide>_<area>_<label>,
# so these carry the real names. YF2's is the case that matters: its fastq prefix
# says B1 and the instrument says the section sat in D1, because in the prefix
# that field is a SAMPLE NAME. Taken from the real SC2630 delivery.
IMG="${HERE}/visium_images"
rm -rf "$IMG"
mkdir -p "$IMG"
: > "${IMG}/CAVG10243R1_2026-10-05_15-14-00_2026-10-05_14-33-36_H1-DYCH77W_A1_A2.tif"
: > "${IMG}/CAVG10243R1_2026-10-05_15-14-00_2026-10-05_14-33-36_H1-DYCH77W_D1_B1.tif"
: > "${IMG}/YF1_microscope.jpg"
printf 'probe_id,gene_id\nENSCAFP0001,ENSCAFG0001\n' > "${IMG}/probe_set.csv"
printf '{}\n' > "${IMG}/YF1_alignment.json"

# YF1's slide is in the sheet, YF2's is not: --slide covers a run that sat on
# one slide, which is the usual case.
cat > "${HERE}/sample-sheet-spatial.csv" <<EOF
label,fastqs,slide,area,cytaimage,image,probeset,alignment
YF1,${YF1a},H1-DYCH77W,A1,CAVG10243R1_2026-10-05_15-14-00_2026-10-05_14-33-36_H1-DYCH77W_A1_A2.tif,YF1_microscope.jpg,probe_set.csv,YF1_alignment.json
YF1,${YF1b},H1-DYCH77W,A1,CAVG10243R1_2026-10-05_15-14-00_2026-10-05_14-33-36_H1-DYCH77W_A1_A2.tif,YF1_microscope.jpg,probe_set.csv,YF1_alignment.json
YF2,${YF2},,D1,CAVG10243R1_2026-10-05_15-14-00_2026-10-05_14-33-36_H1-DYCH77W_D1_B1.tif,,,
EOF

# no image column filled in at all: spaceranger has nothing to align to, so the
# sheet must be rejected before anything is submitted
cat > "${HERE}/sample-sheet-spatial-noimage.csv" <<EOF
label,fastqs,slide,area
YF2,${YF2},H1-DYCH77W,B1
EOF

# a slide serial with no capture area
cat > "${HERE}/sample-sheet-spatial-noarea.csv" <<EOF
label,fastqs,slide,cytaimage
YF2,${YF2},H1-DYCH77W,CAVG10243R1_2026-10-05_15-14-00_2026-10-05_14-33-36_H1-DYCH77W_D1_B1.tif
EOF

# the serial was lost, so one of spaceranger's keywords is used instead and the
# capture area has to be left out
cat > "${HERE}/sample-sheet-spatial-unknown.csv" <<EOF
label,fastqs,slide,cytaimage
YF2,${YF2},visium-2,CAVG10243R1_2026-10-05_15-14-00_2026-10-05_14-33-36_H1-DYCH77W_D1_B1.tif
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
