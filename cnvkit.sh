#!/bin/bash

#set -eo pipefail

## Variables
# Date
NOW=$(date +"%d\%m\%Y")

# Directories
BAM="/home/idengene/Exoma/aln/bam_files"
ALL_BED="/home/idengene/CNV/GenRef/banco"
BED="/home/idengene/CNV/results/cnvkit/bed_files"
OUT_DIR="/home/idengene/CNV/results/cnvkit/cnr_files"
ANNOT_FILES="/home/idengene/CNV/results/cnvkit/annot_files"

# Files
R_SCRIPT="/home/idengene/scripts/./cnvkit.R"
BANCO="python /home/idengene/scripts/banco.py"
EXOME_DEPTH="/home/idengene/scripts/./ExomeDepth.R"
TESTES="/home/idengene/Exoma/aln/bam_files/*.Tumor.bam"
NORMAL="/home/idengene/Exoma/aln/bam_files/*.Normal.bam"
CNVKIT="python2 /home/idengene/CNV/bin/cnvkit-0.9.6/cnvkit.py"

# References
BANCO_N="/home/idengene/CNV/GenRef/banco.bed"
CNV_REF="/home/idengene/CNV/GenRef/ref_CNV.cnn"
REFFLAT="/home/idengene/CNV/GenRef/refFlat.txt"
ACCESS="/home/idengene/CNV/GenRef/my_access.bed"
TARGETS="/home/idengene/CNV/GenRef/targets_105V2.bed"
GEN_REF="/home/idengene/Exoma/GenRef/hg19/chromFa/hg19.fa"

################################################################################## PREPARE BED FILE #######################################################################

# cnvkit.py target Capture.bed --split -a 12 -o my_targets.bed
# cnvkit.py antitarget my_targets.bed -a 1500 -g data/access-5k-mappable.hg19.bed -o my_antitargets.bed

############################################################################## USE CNVKIT ##################################################################################

rm $CNV_REF &&
# rm $BED/*.bed &&
rm $OUT_DIR/*.* &&
rm $ANNOT_FILES/*.* &&

## This command will build all references
$CNVKIT batch --normal $NORMAL --targets $TARGETS --annotate $REFFLAT --fasta $GEN_REF --access $ACCESS --output-reference $CNV_REF --output-dir $OUT_DIR &&

rename 'Normal' Tumor $BAM/*.ba* &&

## Exclude Antitargets regions from ref file
sed -i '/Antitarget/d' $CNV_REF &&

## Run analysis on new Samples using the same references
$CNVKIT batch $TESTES -r $CNV_REF -d $OUT_DIR &&

## Call cn's

cd $OUT_DIR &&

for i in $(ls *.Tumor.cnr | rev | cut -c 11- | rev | uniq); 
    do 
        $CNVKIT call ${i}.Tumor.cnr -y -m clonal -o ${i}.call.cnr && 
        sed -i '/Antitarget/d' ${i}.call.cnr;
done;

$R_SCRIPT &&

cd $ALL_BED &&

# Exclude '"' and ".call.cnr" from R output
for i in $(ls *.bed | rev | cut -c 5- | rev | uniq); 
    do
        sed -i 's/"//g' ${i}.bed &&
        sed -i 's/chr//g' ${i}.bed;
done;

$BANCO &&

rm $ALL_BED/all.amp.bed &&

cd $BED &&

rename '.call.cnr.bed' .bed *.bed &&

# Exclude '"' and ".call.cnr" from R output
for i in $(ls *.bed | rev | cut -c 5- | rev | uniq); 
    do
        sed -i 's/"//g' ${i}.bed &&
        sed -i 's/.call.cnr//g' ${i}.bed;
done;

## Match with Database

for i in $(ls *.bed | rev | cut -c 5- | rev | uniq); 
    do
        awk 'FNR==NR{a[$2,$5]=$4;next} ($2,$4) in a{print $0,a[$2,$4]}' $BANCO_N ${i}.bed > ${i}.final.bed &&
        rm ${i}.bed &&        
#        awk 'FNR==NR{a[$2]=$7;next}{print $0,a[$2]?a[$2]:"0"}' $BANCO_N ${i}.final.bed > ${i}.bed &&
#        rm ${i}.final.bed &&
        sed -i 's/ /\t/g' ${i}.final.bed
done;

## PERFORM ANNOTATIONS FOR .bed FILES
set $ANNOTSV="/home/idengene/CNV/bin/AnnotSV_2.2/bin" &&
export ANNOTSV=/home/idengene/CNV/bin/AnnotSV_2.2 &&

for i in $(ls *.final.bed | rev | cut -c 11- | rev | uniq); 
    do
        $ANNOTSV/bin/AnnotSV/AnnotSV.tcl -SVinputFile ${i}.final.bed -SVinputInfo 1 -typeOfAnnotation split -svtBEDcol 4 -genomeBuild GRCh37 -metrics fr -outputDir $ANNOT_FILES -outputFile ${i}.tsv;

done;

############################################################################### FIM #######################################################################################

