#!/bin/bash

set -eo pipefail

############################################################################# VARIABLES ############################################################################

# OUT FILES
CAROL_OUT="/home/idengene/carol/VEPs"
BAM="/home/idengene/Exoma/aln/bam_files"
VEP_OUT="/home/idengene/Exoma/results/VEP"
GATK_OUT="/home/idengene/Exoma/results/GATK"
TAPES_OUT="/home/idengene/test_tools/tapes/VCFs_teste"

# REFERENCES
VEP_CACHE="/home/idengene/.vep"
VEP_PLUGINS="/home/idengene/.vep/Plugins/"
DBSNP="/home/idengene/Exoma/GenRef/All.vcf.gz"
CLINVAR="/home/idengene/Exoma/GenRef/clinvar.vcf.gz"
REF_GATK="/home/idengene/Exoma/GenRef/hg19/chromFa/hg19.fa"

# TOOLS
VEP="/home/idengene/Exoma/bin/ensembl-vep/./vep"
TAPES="python3 /home/idengene/test_tools/tapes/tapes.py"
FILTER_VEP="/home/idengene/Exoma/bin/ensembl-vep/./filter_vep"

########################################################################### START PIPELINE #########################################################################

echo "Cleaning output folders"

rm $VEP_OUT/* &
rm $TAPES_OUT/*.vcf &
rm -rf $TAPES_OUT/out_files/*/ &&
export BCFTOOLS_PLUGINS=/home/idengene/Exoma/bin/bcftools/plugins/ &&

cd $GATK_OUT &&

for i in $(ls *.vcf.gz | rev | cut -c 8- | rev | uniq);
    do
	echo
	echo "Annotating ${i} for tapes"
	echo

    ## Annotate files for Tapes
	$VEP -i ${i}.vcf.gz --assembly GRCh37 --format vcf -o $TAPES_OUT/${i}.vep.vcf --vcf --sift p --polyphen p --ccds --hgvs --symbol --numbers --regulatory --canonical --protein --biotype --gene_phenotype --af_gnomad --max_af --pubmed --variant_class --hgvs --fasta $REF_GATK --dir $VEP_CACHE --no_stats --use_given_ref --force_overwrite --offline --cache --merged --flag_pick --fork 20 --verbose --plugin MaxEntScan,/home/idengene/Exoma/GenRef/fordownload,SWA,NCSS --plugin SingleLetterAA --plugin dbscSNV,/home/idengene/Exoma/GenRef/dbscSNV1.1_GRCh37.txt.gz --plugin Carol --plugin dbNSFP,/home/idengene/Exoma/GenRef/dbNSFPv3.5a_hg19.gz,gnomAD_genomes_AF,gnomAD_exomes_AF,CADD_phred,FATHMM_converted_rankscore,clinvar_clnsig,clinvar_golden_stars,Interpro_domain,SIFT_score,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_score,MetaSVM_pred,MetaLR_pred,M-CAP_pred,fathmm-MKL_coding_pred,GenoCanyon_score,GERP++_RS --plugin gnomADc,/home/idengene/Exoma/GenRef/gnomAD --plugin gnomADc,/home/idengene/Exoma/GenRef/gnomADg --custom $CLINVAR,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN --custom $DBSNP,dbSNP,vcf,exact,0,ID &

	echo
	echo "Creating VEP files for ${i}"
	echo

    ## Annotate files for Carol
        $VEP -i ${i}.vcf.gz --assembly GRCh37 --format vcf -o $VEP_OUT/${i}.vep --tab --fields "ClinVar_CLNSIG,gnomAD_AF,Consequence,SYMBOL,HGVSc,HGVSp,Existing_variation,Location,${i}_POS,${i}_REF,${i}_ALT,${i}_DP,${i}_SOR,Codons,HGVSg,REVEL,CLIN_SIG,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,MAX_AF_POPS,SIFT,PolyPhen,CAROL,MES-NCSS_downstream_donor,MES-NCSS_upstream_acceptor,MES-NCSS_upstream_donor,MES-SWA_acceptor_alt,MES-SWA_acceptor_diff,MES-SWA_acceptor_ref,MES-SWA_acceptor_ref_comp,MES-SWA_donor_alt,MES-SWA_donor_diff,MES-SWA_donor_ref,MES-SWA_donor_ref_comp,MaxEntScan_alt,MaxEntScan_diff,MaxEntScan_ref,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_rankscore,BayesDel_noAF_pred,CADD_phred,VEST4_score,ClinVar_CLNREVSTAT,ClinVar_CLNDN,ClinVar,Feature,Feature_type,cDNA_position,CDS_position,Protein_position,Amino_acids,DISTANCE,IMPACT,STRAND,VARIANT_CLASS,BIOTYPE,CANONICAL,EXON,INTRON,SOMATIC,Allele,PUBMED" --sift p --polyphen p --ccds --hgvs --hgvsg --symbol --keep_csq --numbers --regulatory --canonical --protein --biotype --gene_phenotype --af_gnomad --max_af --pubmed --variant_class --fasta $REF_GATK --cache --dir $VEP_CACHE --dir_plugin $VEP_PLUGINS --no_stats --use_given_ref --force_overwrite --offline --merged --flag_pick --fork 20 --verbose --plugin MaxEntScan,/home/idengene/Exoma/GenRef/fordownload,SWA,NCSS --plugin dbscSNV,/home/idengene/Exoma/GenRef/dbscSNV1.1_GRCh37.txt.gz --plugin Carol --plugin SubsetVCF,file=${i}.vcf.gz,name=${i},filter=0,fields=DP%SOR --plugin dbNSFP,/home/idengene/Exoma/GenRef/dbNSFP4.1a_grch37.gz,gnomAD_genomes_AF,gnomAD_exomes_AF,CADD_phred,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_rankscore,BayesDel_noAF_pred,FATHMM_converted_rankscore,clinvar_clnsig,Interpro_domain,SIFT_score,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_score,MetaSVM_pred,MetaLR_pred,M-CAP_pred,fathmm-MKL_coding_pred,GenoCanyon_score,GERP++_RS,VEST4_score --plugin gnomADc,/home/idengene/Exoma/GenRef/gnomAD --plugin REVEL,/home/idengene/Exoma/GenRef/new_tabbed_revel.tsv.gz --custom $CLINVAR,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN --custom $DBSNP,dbSNP,vcf,exact,0,ID &&

	echo
	echo "Filtering ${i}.vep and including QUAL and FILTER columns"
	echo

    ## Keep only canonical isoforms
	$FILTER_VEP -i $VEP_OUT/${i}.vep --format tab --filter "Feature in /home/idengene/Exoma/GenRef/motifs_list.txt" --force_overwrite -o $VEP_OUT/${i}.mod.vep &&
	rm $VEP_OUT/${i}.vep &&
    	
    # Inclui colunas QUAL, FILTER e VAF no arquivo final.vep
	bcftools +fill-tags ${i}.vcf.gz -Oz -o ${i}_MAF.vcf.gz -- -t VAF &&
	gunzip ${i}_MAF.vcf.gz &&
	sed -i '/^##/d' ${i}_MAF.vcf &&
	sed -i '/^##/d' $VEP_OUT/${i}.mod.vep &&
	vcftools --vcf ${i}_MAF.vcf --extract-FORMAT-info VAF --out ${i} &&
	awk '{ if ($65 == "insertion" || $65 == "deletion") $9 = $9-1 ; print}' $VEP_OUT/${i}.mod.vep > $VEP_OUT/tmp && 
	mv $VEP_OUT/tmp $VEP_OUT/${i}.mod.vep &&
	awk -v OFS="\t" 'NR==FNR{a[$2]=$7;next}{$6=$6 "\t"(a[$9]?a[$9]:"FILTER")}1' ${i}_MAF.vcf $VEP_OUT/${i}.mod.vep > $VEP_OUT/${i}.final.vep &&
	awk -v OFS="\t" 'NR==FNR{a[$2]=$6;next}{$6=$6 "\t"(a[$10]?a[$10]:"QUAL")}1' ${i}_MAF.vcf $VEP_OUT/${i}.final.vep > $VEP_OUT/${i}.mod.vep &&
	awk -v OFS="\t" 'NR==FNR{a[$2]=$3;next}{$14=$14 "\t"(a[$11]?a[$11]:"VAF")}1' ${i}.VAF.FORMAT $VEP_OUT/${i}.mod.vep > $VEP_OUT/${i}.final.vep &&	
	rm ${i}.log &&
	rm ${i}_MAF.vcf &&
	rm ${i}.VAF.FORMAT &&
	rm $VEP_OUT/${i}.mod.vep
done;

rm $CAROL_OUT/*.* &&

grep -H -w -e "rs17107315" -e "pathogenic" $VEP_OUT/*.final.vep > $CAROL_OUT/pathogenic_variants.txt &&

cd $TAPES_OUT &&

    ## Run Tapes on annotated files
for i in $(ls *.vep.vcf | rev | cut -c 9- | rev | uniq); 
    do 
	awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${i}.vep.vcf > ${i}.chr.vcf && rm ${i}.vep.vcf &&
        $TAPES sort -i $TAPES_OUT/${i}.chr.vcf -o $TAPES_OUT/out_files/${i}/ --tab --by_sample --by_gene -a hg19 
done;

cp $VEP_OUT/* $CAROL_OUT &&
cp $TAPES_OUT/out_files/*/*.txt $CAROL_OUT &&
cp $GATK_OUT/*.gz $CAROL_OUT &&
gunzip $CAROL_OUT/*.gz &&
cp $BAM/gender.txt $CAROL_OUT &&

echo
echo "Pipeline Finished"
echo

########################################################################### END PIPELINE #########################################################################

