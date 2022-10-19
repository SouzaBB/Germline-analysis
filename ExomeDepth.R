#!/usr/bin/Rscript

#################################################################### Use ExomeDepth ####################################################################

## Load library
library(ExomeDepth)
library(GenomicRanges) 

## get the annotation datasets to be used later
data(exons.hg19) #included
exons.hg19.GRanges <- GRanges(seqnames = exons.hg19$chromosome,IRanges(start=exons.hg19$start,end=exons.hg19$end),names = exons.hg19$name) 

## Load all BAM files
path <- "/home/idengene/Exoma/aln/bam_files/"
files = list.files(path = path, pattern = "*.Tumor.bam", full.names=FALSE)

## Modify targets file to work with ExomeDepth
segments <- read.table("/home/idengene/CNV/GenRef/targets_105V2.bed", sep="\t", as.is=TRUE)
colnames(segments) <-c("chromosome","start","end", "gene")
as.character(segments$chromosome) -> segments$chromosome

## Create count data from BAM files
my.counts <- getBamCounts(bed.frame = segments, bam.files = files, include.chr = FALSE, referenceFasta = '/home/idengene/Exoma/GenRef/hg19/chromFa/hg19.fa')

ExomeCount.dafr <- as(my.counts, 'data.frame')
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome), pattern = 'chr', replacement = '')  ##remove the annoying chr letters
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*bam')])
nsamples <- ncol(ExomeCount.mat)

as.character(seq(1:length(ExomeCount.dafr$end))) -> exons

for(i in 1:nsamples) {
  my.choice <- select.reference.set(test.counts = ExomeCount.mat[,i], reference.counts = ExomeCount.mat[,-i],bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,n.bins.reduced = 10000)
  print(my.choice)
  my.reference.selected <- apply(X = ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE], MAR = 1,FUN = sum)
  all.exons <- new('ExomeDepth',test = ExomeCount.mat[,i],reference = my.reference.selected,formula = 'cbind(test, reference) ~ 1')
  all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4, chromosome = ExomeCount.dafr$chromosome, start = ExomeCount.dafr$start, end = ExomeCount.dafr$end, name=exons)
  all.exons <- AnnotateExtra(x = all.exons,reference.annotation = exons.hg19.GRanges,min.overlap = 0.0001,column.name = 'exons.hg19') 
  output.file <- paste(colnames(ExomeCount.mat)[i], '_ExomeDepth.bed', sep = '')
  write.table(file = output.file, x = all.exons@CNV.calls, row.names = FALSE, sep="\t", dec=".")
}

######################################################## Convert ExomeDepth output to .BED file ########################################################

path <- "/home/idengene/Exoma/aln/bam_files/"
files = list.files(path = path, pattern = "*.bed", full.names=FALSE)
for(file in files) {
	perpos <- which(strsplit(file, "")[[1]]==".")
	assign(
	gsub(" ","",substr(file, 1, perpos-1)), 
	read.csv(paste(path,file,sep=""), header=T, sep="\t"))

}

out_BED = function(file) {
	file_df <- read.csv(file, header=T, sep="\t")

	# Calculate log2 ratio and absolute cn
    	file_df$log2 <- round(log2(file_df$reads.ratio), digits=3)
    	file_df$cn <- round(2*(2^(file_df$log2)))	

	# Re-order columns for better view
	ordered <- file_df[c(7,5,6,3,15,13,9,11,10,12)]

	# Convert "ordered" to .bed file
	names(ordered)[1] <- "#Chrom"
	names(ordered)[2] <- "Start"
	names(ordered)[3] <- "End"
	names(ordered)[4] <- "SV_type"
    	names(ordered)[6] <- "Gene_exon"	
    	names(ordered)[8] <- "reads.Sample"
    	names(ordered)[9] <- "reads.Reference"
	names(ordered)[10] <- "reads.ratio" 	

 	# Rename 'deletion' and 'duplication' to 'DEL' and 'AMP'    
    	ordered$SV_type <- ifelse(ordered$SV_type == "deletion", "DEL", "AMP")
	ordered$BF <- as.numeric(ordered$BF)
	ordered$Prediction <- ifelse(ordered$BF >= 175, "Suspect", "Probably Benign")

	# Save "ordered" as .bed file
	output_file = sprintf("/home/idengene/CNV/results/cnvkit/bed_files/%s", file)
	write.table(ordered, file=output_file, row.names=F, col.names=T, sep="\t", dec=",")
	#output_file = sprintf("/home/idengene/carol/CNVs/%s", file)
        #write.table(ordered, file=output_file, row.names=F, col.names=T, sep="\t", dec=",")

	# Create VCF format columns
	#file_vcf <- ordered	
	#file_vcf$ID <- "." 
	#file_vcf$REF <- "."
	#file_vcf$QUAL <- "."
	#file_vcf$FILTER <- "."
	#file_vcf$INFO <- paste0("SVTYPE=",file_vcf$SV_type,";END=",file_vcf$End)
	#file_vcf$FORMAT <- "."
	#final_vcf <- file_vcf[c(1,2,12,13,4,14,15,16,17)]
	
	# Rename columns
	#names(final_vcf)[1] <- "CHROM"
	#names(final_vcf)[2] <- "POS"
	#names(final_vcf)[5] <- "ALT"

	# Export VCF format
	#output_file = sprintf("/home/idengene/CNV/results/cnvkit/vcf_files/%s", file)
	#write.table(final_vcf, file=output_file, row.names=F, col.names=T, sep="\t")

}
	print(files)
	for (file in files) {
	    out_BED(file)
}

########################################################################################################################################################

