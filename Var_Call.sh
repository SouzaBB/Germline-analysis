#!/bin/bash

#set -eo pipefail

############################################################################# VARIABLES ############################################################################

# ALIGNMENTS 
BAM="/home/idengene/Exoma/aln/bam_files"

# OUT FILES
GATK_OUT="/home/idengene/Exoma/results/GATK"
FINAL_VCF="/home/idengene/Exoma/results/Final_VCF"

# REFERENCES
DBSNP="/home/idengene/Exoma/GenRef/All.vcf.gz"
TARGETS="/home/idengene/CNV/GenRef/targets_105V2.bed"
DICT="/home/idengene/Exoma/GenRef/hg19/chromFa/hg19.dict"
REF_GATK="/home/idengene/Exoma/GenRef/hg19/chromFa/hg19.fa"

# TOOLS
GATK="java -jar -Xmx20G /home/idengene/Exoma/bin/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar"

########################################################################### START PIPELINE #########################################################################

cd $BAM &&
rm $GATK_OUT/*.* &&

## Variant call
for i in $(ls *.Tumor.bam | rev | cut -c 11- | rev | uniq); 
    do
        
    ## Variant call GATK
        $GATK HaplotypeCaller -R $REF_GATK -I ${i}.Tumor.bam -O $GATK_OUT/${i}.1.vcf --sequence-dictionary $DICT --dbsnp $DBSNP -stand-call-conf 30 --min-pruning 3 --add-output-vcf-command-line false -L $TARGETS &&

	$GATK VariantFiltration -V $GATK_OUT/${i}.1.vcf -filter "QD < 2.0" --filter-name "LowQD" -filter "QUAL < 30.0" --filter-name "LowQual" -filter "FS > 60.0" --filter-name "SB" -filter "MQ < 20.0" --filter-name "LowMQ" -filter "DP < 100.0" --filter-name "LowDP" --add-output-vcf-command-line false -O $GATK_OUT/${i}.2.vcf &&

	bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' $GATK_OUT/${i}.2.vcf -o $GATK_OUT/${i}.vcf &&

	bgzip -f $GATK_OUT/${i}.vcf > $GATK_OUT/${i}.vcf.gz &&

	tabix -p vcf $GATK_OUT/${i}.vcf.gz &&

	rm $GATK_OUT/${i}.1.vcf* &&

	rm $GATK_OUT/${i}.2.vcf*

done;

## Get Coverage on Y chromosome
for i in $(ls *.Tumor.bam | rev | cut -c 11- | rev | uniq); 
       do 
               samtools depth -a -r Y:2654985-2655305 ${i}.Tumor.bam | awk '{sum+=$3} END { print "Average of sample '${i}' = ", sum/NR }' > ${i}_gender.txt
done;

more *_gender.txt > gender.txt
