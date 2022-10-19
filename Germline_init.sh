#!/bin/bash

############################################################################################################################################################

BAM="/home/idengene/Exoma/aln/bam_files"
BASESPACE="/home/idengene/basespace/./bs"
CNV="/home/idengene/scripts/CNV_Germline.sh"
READS_FOLDER="/home/idengene/basespace/Samples"
VCF_FOLDER="/home/idengene/basespace/Samples/VCF"
ANOT="bash /home/idengene/scripts/Anotação_Carol.sh"

############################################################################################################################################################

$BASESPACE project list

echo "Is your project on the list?"
read -r -p "'Yes' or 'No' [Y/N] = " response

if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]  
then
	read -r -p "Give me the extension you want [fastq|VCF]: " extension	
	read -r -p "Give me the project name: " project
	if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]] && [[ "$extension" =~ ^([vV][cC][fF]|[VCF])$ ]]
	then
		rm -rf $VCF_FOLDER/*_ds.*/
		rm $VCF_FOLDER/VCFs/*
		$BASESPACE download project --name $project --extension ${extension,,} -o $VCF_FOLDER
		$ANOT			
		echo
	elif [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]] && [[ "$extension" =~ ^([fF][aA][sS][tT][qQ]|[fastq])$ ]]
	then
		rm -rf $READS_FOLDER/*_ds.*/
		rm -rf $READS_FOLDER/reads/*
		rm $BAM/*.ba*
		$BASESPACE download project --name $project --extension ${extension,,}.gz -o $READS_FOLDER 
		$CNV
		echo
	else
		echo "'${extension,,}' is not a valid extension!"
	fi
else
	echo "Need to authenticate on the right account"
	$BASESPACE auth --force
fi

#TODO VCF como minúsculo - DONE
#TODO Execute pipeline from anywhere - DONE
#TODO Identify type of analysis and execute script - DONE
#TODO Add biosample support
#TODO while statement for when user needs to login into another profile

