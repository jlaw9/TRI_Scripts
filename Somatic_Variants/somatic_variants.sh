#!/bin/bash

#04/14/2015 
#Author: JeffL

#Goal: Get the somatic variants from the 3x3 table VCF files

if [ $# -lt 2 ]; then
	echo "USAGE bash get_somatic_variants.sh </path/to/N-1vsT-1> <Sample_Name>"
	exit 8
fi

DIR="$1"
SAMPLE="$2"
mkdir -p "${DIR}/Somatic_Variants"
SV_DIR="${DIR}/Somatic_Variants"
CHR="$3"

#if [ "`find /home/ionadmin/jeff/Lung_Somatic/${SAMPLE}_somatic.vcf 2>/dev/null`" ]; then
#	echo "The somatic variants for $SAMPLE have already been annotated. Skipping"
#	exit 0
#fi

#Somatic Variants scripts DIR
SV_SCRIPTS="/home/ionadmin//TRI_Scripts/Somatic_Variants"

#first extract the somatic variants
grep -E "WT.+HET" ${DIR}/matched_variants${CHR}.csv > ${SV_DIR}/somatic.csv

# copy the variant if there is only one so that Ozlem's scripts will work
if [ "`head ${SV_DIR}/somatic.csv | wc -l `" == "1" ]; then
	grep -E "WT.+HET" ${DIR}/matched_variants${CHR}.csv >> ${SV_DIR}/somatic.csv
fi


#now get the variant info for the somatic variants by matching somatic.csv and VCF2${CHR}_Final.vcf, then outputting somatic.vcf
python2.7 ${SV_SCRIPTS}/getVCF_v2.py ${SV_DIR}/somatic.csv ${DIR}/VCF2${CHR}_Final.vcf ${SV_DIR}/somatic.vcf

# Annotate the somatic variants with table annovar. 
# I know we updated some of the annovar code in the file /home/ionadmin//TRI_Scripts/annovar/annotate_variation.pl starting on line 1552
/home/ionadmin//TRI_Scripts/annovar/convert2annovar.pl --format vcf4old --snpqual 0 ${SV_DIR}/somatic.vcf > ${SV_DIR}/somatic_input.vcf
/home/ionadmin//TRI_Scripts/annovar/table_annovar.pl \
	${SV_DIR}/somatic_input.vcf \
	/rawdata/software/annovar/humandb_ucsc/ \
	--outfile ${SV_DIR}/annovar_somatic_table \
	--buildver hg19 \
	--protocol refGene,1000g2012apr_all,snp129NonFlagged,snp137NonFlagged,cosmic68,ljb23_all \
	--operation g,f,f,f,f,f \
	--remove \
	--nastring . 


# Add the allele frequencies and such to the annotated somatic variants

python2.7 ${SV_SCRIPTS}/updateAnnovarTable_v3.py \
	${SV_DIR}/somatic.csv \
	${SV_DIR}/annovar_somatic_table.hg19_multianno.txt \
	${SV_SCRIPTS}/transcriptList.txt \
	${SV_DIR}/annovar_summary.txt 

# Ozlem's notes about the transcriptList.txt file: I believe this is the canonical transcript list for whole exome (I wish I had time to double-check though!) so please keep in mind that I did not have a chance to double-check: /results/ozlem/pnet_merge_analysis/transcriptList.txt.        
# Just noticed that somatic.vcf has duplicate entries, will need to correct this 

# copy the final file back to my dir
echo "copying ${SV_DIR}/annovar_summary.txt to /home/ionadmin/jeff/Lung_Somatic/"
mkdir -p /home/ionadmin/jeff/Lung_Somatic/
cp ${SV_DIR}/annovar_summary.txt /home/ionadmin/jeff/Lung_Somatic/${SAMPLE}_somatic.vcf
if [ $? != 0 ]; then
	echo "Failed!"
	exit 1
fi
echo "Finished"


