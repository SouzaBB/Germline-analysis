#!/bin/bash

#set -eo pipefail

############################################################################# VARIABLES ############################################################################

# ALIGNMENTS 
BAM="/home/idengene/Exoma/aln/bam_files"
LOG="/home/idengene/Exoma/aln/log_files"
SAM="/home/idengene/Exoma/aln/sam_files"

# CNV call scripts
CNV="/home/idengene/scripts/./cnvkit.sh"
EXOME_DEPTH="Rscript /home/idengene/scripts/ExomeDepth.R"

# INPUT FILES
SAMPLES="/home/idengene/samples/"
FASTQS="/home/idengene/basespace/Samples"
READS="/home/idengene/basespace/Samples/reads/"

# OUT FILES
ALU_OUT="/home/idengene/CNV/results/MELT"
FASTQC_OUT="/home/idengene/Exoma/results/FastQC"
CNV_OUT="/home/idengene/CNV/results/cnvkit/bed_files"

# REFERENCES
DBSNP="/home/idengene/Exoma/GenRef/All.vcf.gz"
REF="/home/idengene/Exoma/GenRef/hg19/chromFa/hg19"
TARGETS="/home/idengene/CNV/GenRef/targets_105V2.bed"
DICT="/home/idengene/Exoma/GenRef/hg19/chromFa/hg19.dict"
REF_GATK="/home/idengene/Exoma/GenRef/hg19/chromFa/hg19.fa"
ALU_MELT="/home/idengene/CNV/bin/MELTv2.2.0/me_refs/1KGP_Hg19/ALU_MELT.zip"
GENES_BED="/home/idengene/CNV/bin/MELTv2.2.0/add_bed_files/1KGP_Hg19/hg19.genes.bed"

# TOOLS
BWA="/home/idengene/Exoma/bin/bwa/./bwa"
FASTQC="/home/idengene/Exoma/bin/FastQC/./fastqc"
MELT="java -jar -Xmx15G /home/idengene/CNV/bin/MELTv2.2.0/MELT.jar"
PICARD="java -jar -Xmx15G /home/idengene/Exoma/bin/picard-2.23.0/picard.jar"
GATK="java -jar -Xmx15G /home/idengene/Exoma/bin/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar"

# Variant Calling and Annotation
VAR_CALL="bash /home/idengene/scripts/Var_Call.sh"
VAR_ANNOT="bash /home/idengene/scripts/Var_annot.sh"

######################################################################## CREATE REFERENCES #########################################################################
########################################################################### DO ONLY ONCE ###########################################################################

## echo "Creating references for BWA e PICARD"

# $BWA index -a bwtsw -p hg19 $REF &&
# samtools faidx $REF_GATK &&
# $PICARD CreateSequenceDictionary -R $REF_GATK -O $DICT

## echo "All references created"
## echo

######################################################################## CREATE DIRECTORIES ########################################################################

echo "Creating output directories"

mkdir -p /home/idengene/Exoma/aln/{bam_files,sam_files} &&
mkdir -p /home/idengene/Exoma/results/{FastQC,Final_VCF,GATK,VEP} &&

echo "All directories created"
echo

########################################################################### START PIPELINE #########################################################################

echo "Preparing files"
echo

mv $FASTQS/*/*.fastq.gz $READS &&

cd $READS &&
fastq_lane_merging.sh &&
mv *_ME_* $SAMPLES &&

rm $CNV_OUT/*.* &&

echo
echo "Output folder cleaned"
echo

cd $SAMPLES &&

echo
echo "Starting Pipeline in $(ls -1 *_R1* | wc -l) Samples"
echo

## Rename all fastq files
for i in $(ls *.fastq.gz | cut -d'-' -f1 | uniq); 
    do 
        mv ${i}-*_R1_*.fastq.gz ${i}_R1.fastq.gz && 
        mv ${i}-*_R2_*.fastq.gz ${i}_R2.fastq.gz
done;

## Alignment, SAM to BAM conversion and mark duplicates
for i in $(ls *.fastq.gz | rev | cut -c 13- | rev | uniq); 
    do

    ## QC 	
#        $FASTQC -t 14 ${i}_R1.fastq.gz ${i}_R2.fastq.gz -o $FASTQC_OUT/ &&

    ## alignment
        $BWA mem -M -t 14 -R "@RG\tID:idengene\tLB:AGILENT\tSM:${i}\tPL:ILLUMINA" $REF ${i}_R1.fastq.gz ${i}_R2.fastq.gz > $SAM/${i}.aln.sam &&
        echo
            
    ## Clean Sam
        $PICARD CleanSam I=$SAM/${i}.aln.sam O=$SAM/${i}.sam &&
        rm $SAM/${i}.aln.sam &&
        echo
        
    ## SAM to BAM conversion
        $PICARD SortSam I=$SAM/${i}.sam O=$BAM/${i}.bam SO=coordinate CREATE_INDEX=true &&
        rm $SAM/${i}.sam &&
        echo

    ## mark duplicates
        $PICARD MarkDuplicates I=$BAM/${i}.bam O=$BAM/${i}.Normal_1.bam METRICS_FILE=$BAM/${i}.bam.metrics CREATE_INDEX=true &&
        rm $BAM/${i}.{bam,bai} &&
        rm ${i}_R{1,2}.fastq.gz 
        echo
      
    ## Base Recalibration
        $GATK BaseRecalibrator -R $REF_GATK -I $BAM/${i}.Normal_1.bam --known-sites $DBSNP -O $BAM/${i}.recal.data.csv -L $TARGETS &&
        $GATK ApplyBQSR -R $REF_GATK -I $BAM/${i}.Normal_1.bam -bqsr $BAM/${i}.recal.data.csv -O $BAM/${i}.Normal.bam
	rm $BAM/${i}.Normal_1.{bam,bai}
	echo  
done;

cd $BAM &&

## Run CNV analysis
$CNV && 
$EXOME_DEPTH && 
rm $BAM/X* &&
$VAR_CALL && 
$VAR_ANNOT &&

# Rename bam files
rename 'Tumor.bai' Tumor.bam.bai *.bai &&	

## Alu call  
for i in $(ls *.Tumor.bam | rev | cut -c 11- | rev | uniq); 
    do
   	mkdir -p $ALU_OUT/${i}/
	$MELT Single -h $REF_GATK -ac -bamfile $BAM/${i}.Tumor.bam -n $GENES_BED -t $ALU_MELT -w $ALU_OUT/${i} -c 35
	rm $BAM/*.disc*
	rm $BAM/*.fq
done;

################################################################################# END ##############################################################################

