#!/usr/bin/Rscript

# Clear workspace
rm(list=(ls()))

path <- "/home/idengene/CNV/results/cnvkit/cnr_files/"
files = list.files(path = path, pattern = "*.call.cnr", full.names=FALSE)
for(file in files) {
	perpos <- which(strsplit(file, "")[[1]]==".")
	assign(
	gsub(" ","",substr(file, 1, perpos-1)), 
	read.csv(paste(path,file,sep=""), header=T, sep="\t"))

}

s_ref <- read.csv("/home/idengene/CNV/GenRef/ref_CNV.cnn", header=T, sep="\t")
#s_ref["depth.ref.norm"] <- round(s_ref["depth"]/mean(s_ref[["depth"]]), digits=2)

mod_CNV = function(file) {
	file_df <- read.csv(file, header=T, sep="\t")

	# Normalize $depth by mean
   	#file_df[sprintf("depth.%s.norm", file)] <- round(file_df[["depth"]]/mean(file_df[["depth"]]), digits=2)

	# Merge both files by "start" position
	merged <- merge(file_df, s_ref, by="start", suffixes=c(sprintf(".%s", file), ".ref"), all=TRUE)

	# Round "log2" column
	log2_col_name = sprintf("log2.%s", file)
	merged[log2_col_name] <- round(merged[[log2_col_name]], digits=1)

	# re-calculate "cn" based on log2 correction
	merged["cn"] <- round(2*(2^(merged[[log2_col_name]])))

	# Subset file with all "cn" values that are not 2
	alt_cn <- subset(merged, merged[["cn"]] != 2)

	# Create new data with columns of interest
	alt_cns <- as.data.frame(alt_cn[, c(1:4,6:8,13)])

	# Re-order columns for better view
	alt_cns <- alt_cns[c(2,1,3,4,5,7,6,8)]

	# Calculate ratio between coverages
        #alt_cns["depth.ratio.norm"] <- round(alt_cns[[sprintf("depth.%s.norm", file)]] / alt_cns[["depth.ref.norm"]], digits=2)
	alt_cns["depth.ratio"] <- round(alt_cns[[sprintf("depth.%s", file)]] / alt_cns[["depth.ref"]], digits=2)

	## Function to call for DUP or DEL.  
	alt_cns["SV_type"] <- ifelse(alt_cns$cn < 2, "DEL", "AMP")

    	# Convert "alt.cns" to .bed file
	full <- alt_cns[c(1,2,3,10,5,4,6,7,8,9)]
	names(full)[1] <- "#Chrom"
	names(full)[2] <- "Start"
	names(full)[3] <- "End"
	names(full)[4] <- "SV_type"
	names(full)[6] <- "gene"
    
    	full["weight"] <- round(full[["weight"]], digits = 2)
    	full <- full[order(full$"#Chrom"),]

	# Save "full" as .bed file
	output_file = sprintf("/home/idengene/CNV/results/cnvkit/bed_files/%s.bed", file)
	write.table(full, file=output_file, row.names=F, col.names=T, sep="\t", dec=",")

	# Output files to create database
   	AMP <- subset(full, full[["SV_type"]] == "AMP")
   	DEL <- subset(full, full[["SV_type"]] == "DEL")

    	output_file = sprintf("/home/idengene/CNV/GenRef/banco/%s.amp.bed", file)
	write.table(AMP, file=output_file, row.names=F, col.names=T, sep="\t", dec=",")

    	output_file = sprintf("/home/idengene/CNV/GenRef/banco/%s.del.bed", file)
	write.table(DEL, file=output_file, row.names=F, col.names=T, sep="\t", dec=",")

}
	print(files)
	for (file in files) {
		mod_CNV(file)
}

#####################################################################################################################################################

# pdf("../Rplots/NA14234_APC.pdf", height = 4.0)
# ggplot(NA18949_BRCA1, aes(x=start,y=cn)) + geom_point(aes(color = factor(cn)),size=2) + ylim(0, 4) + geom_hline(alpha=0.5, size=0.8, yintercept=2, color='black') + labs(title="Patient_Name", subtitle="BRCA1", x="Coordenadas", y="Copy_Number", caption = "**Gráfico de resultados gerados por NGS") + scale_colour_manual(name="Número_de_Cópias", labels=c("0" = "Perda", "1" = "Perda", "2" = "Normal", "3" = "Ganho", "4" = "Ganho"), values = c("0" = "red", "1" = "red","2" = "black","3" = "blue", "4" = "blue")) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.caption = element_text(color = "green", face = "italic", size = 8))
# dev.off()

