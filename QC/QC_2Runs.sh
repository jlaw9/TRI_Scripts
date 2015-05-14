#!/bin/bash

# GOAL: Run TVC on both run's bam files to generate a 3x3 QC table. 
# TS version: 4.2
# Author: Jeff L

# QSUB variables

#$ -S /bin/bash
#$ -cwd
#$ -N QC_2Runs
#$ -o qc_sample_out.log
#$ -e qc_sample_out.log
#$ -q all.q
#$ -V


SOFTWARE_ROOT="/rawdata/legos"
QC_SCRIPTS="${SOFTWARE_ROOT}/scripts/QC"
#QC_MASTER_OUT="$QC_SCRIPTS/QC_Out_Files/multiple_runs.csv"
BAM_INDEXER='/opt/picard/picard-tools-current/BuildBamIndex.jar'
VARIANT_CALLER_DIR='/results/plugins/variantCaller'
REF_FASTA='/results/referenceLibrary/tmap-f3/hg19/hg19.fasta'
# Default Options
AMP_COV_CUTOFF=30 # The minimum amount of coverage each amplicon needs to have. Default is 30
GET_CDS_DEPTHS="False" # Variable to get depths of cds as well as the project bed.
CLEANUP="False"

function usage {
cat << EOF
USAGE: bash QC_2Runs.sh 
If the output_dir specified already has the vcf files generated from QC, running TVC and samtools depth will be skipped. 
All options up to --gt_cutoffs are required.
	-h | --help
	-r | --run_dirs <path/to/run1_dir> <path/to/run2_dir> (Normal_Dir should always come before Tumor_Dir) (PTRIM.bam is generated. *amplicon.cov.xls is used, and if the .vcf file is not in the run_dir, this script checks in the tvc_out dir)
	-o | --output_dir <path/to/output_dir> (Output dir will hold all of the generated files. Should follow this pattern: sample1/QC/Run1vsRun2.)
	-jo | --json_out <sample/QC/results_QC.json> (This json file will hold the QC error metrics for this comparison. It should be located in the sample's QC dir and should the pattern *_QC.json so QC_generateSheets.py can find it.)
	-b | --project_bed <path/to/project_bed> 
	-a | --amp_cov_cutoff <Amplicon_Coverage_Cutoff> (Cutoff for # of amplicon reads.) 
	-jp | --json_paras <Json_parameters1> <Json_parameters2> (two json parameters files used by TVC)
	-d | --depth_cutoffs <Depth_Cutoff1> <Depth_Cutoff2> (Variants with depth < cutoff will be filtered out)
	-gt | --gt_cutoffs <WT_Cutoff1> <HOM1_Cutoff1> <WT_Cutoff2> <HOM1_Cutoff2>
	-sb | --subset_bed <path/to/subset_bed> (If the user want to subset out certain genes from the Project_Bed)
	-v | --all_vcfs_dir <path/to/allRun1vsRun2> (In order to save time, skip running TVC and just subset the 718 gene set or chr1 from the all)
	-chr | --subset_chr <chr#> (The chromosome specified here (for example: chr1) will be used to subset the VCF and BAM files)
	-cl | --cleanup (If calling QC_getRunInfo.sh after this script, the PRTIM.bam is needed so DO NOT CALL CLEANUP. Delete temp_files used to create the two Output VCF files, the PTRIM.bam the chr_subset bam files if they were created.)
	-B | --bases <total_eligible_bases> <total_possible_bases> (If these total bases have already been calculated you can include them here)
EOF
exit 8
}
	# option not supported anymore
	#-cd | --get_cds_depths (Normally, samtools depth is run using the subset of the beds specified above. If the cds_bed is specified, and this option is specified, samtools depth will be run twice.)
	#-cb | --cds_bed <path/to/CDS_bed> (This option should only be used if the user wants to run TVC using the Project bed, and then intersect the results with the CDS bed.)
	#-cle | --cleanup_everything (Option not yet implemented: Delete everything but the log and the matched_variants.csv)

# Checks to ensure that the files provided exist and are readable. 
# @param $1 should be a list of files.
function checkFiles {
	files=$1
	for file in ${files[@]}; do
		if [ "$file" != "" ]; then
			if [ ! "`find ${file} -maxdepth 0 2>/dev/null`" ]; then
				echo "-- ERROR -- '${file}' not found or not readable" 1>&2
				exit 4
			fi
		fi
	done
}

# $1: RUN1_DIR, $2: RUN2_DIR
# The RUN1_BAM and RUN2_BAM variabels will be used later.
function setupBAMs {
	# for TS v4.2, PTRIM.bam cannot be used as input to TVC.
	RUN1_BAM=`find ${1}/*.bam -maxdepth 0 -type f | grep -v "PTRIM"`
	RUN2_BAM=`find ${2}/*.bam -maxdepth 0 -type f | grep -v "PTRIM"`
	
	if [ "$CHR" != "" ]; then
		# basename gets only the nave of the progam out of the filepath.
		run1_bam_name=`basename $RUN1_BAM`
		run2_bam_name=`basename $RUN2_BAM`
		# Check to see if the index file exists. If it does, then samtools must have finished already before, and it was indexed.
		# If not, then start by getting on the specified chromosome from samtools
		if [ ! "`find ${1}/${CHR}/${CHR}_${run1_bam_name}.bai -maxdepth 0 2>/dev/null`" ]; then
			rm -rf ${1}/${CHR} 2>/dev/null
			mkdir -p ${1}/${CHR} 2>/dev/null
			echo "	Making RUN1_BAM chr subset: $RUN1_BAM" >> $log
			samtools view -b $RUN1_BAM  "$CHR" > ${1}/${CHR}/${CHR}_${run1_bam_name}
		fi
		if [ ! "`find ${2}/${CHR}/${CHR}_${run2_bam_name}.bai -maxdepth 0 2>/dev/null`" ]; then
			rm -rf ${2}/${CHR} 2>/dev/null
			mkdir -p ${2}/${CHR} 2>/dev/null
			echo "	Making RUN2_BAM chr subset for: $RUN2_BAM" >> $log
			samtools view -b $RUN2_BAM  "$CHR" > ${2}/${CHR}/${CHR}_${run2_bam_name}
		fi
		RUN1_BAM="${1}/${CHR}/${CHR}_${run1_bam_name}"
		RUN2_BAM="${2}/${CHR}/${CHR}_${run2_bam_name}"
		checkBamIndex $RUN1_BAM
		checkBamIndex $RUN2_BAM
	fi
}

# $1: the Run_Dir, $2, the Run_num
function setupVCF {
	if [ "$CHR" != "" ]; then
		if [ "`find ${1}/tvc*_out -maxdepth 0 -type d 2>/dev/null`" ]; then
			checkFiles "${1}/tvc*_out/TSVC_variants.vcf" 
			vcf="${1}/tvc*_out/TSVC_variants.vcf"
		else
			checkFiles "${1}/*.vcf"
			vcf=`find ${1}/*.vcf -type f | head -n 1`
		fi
		# what if there are multiple vcf files in the dir? Just take the first one it finds for now.
		grep "^#" $vcf > ${TEMP_DIR}/Run${2}.vcf
		# use grep to get only the specified chromosome out of the vcf file. -E: regex -P: perl regex
		grep -v "^#" $vcf | grep -P "^${CHR}\t" >> ${TEMP_DIR}/Run${2}.vcf
		bgzip -c ${TEMP_DIR}/Run${2}.vcf > ${TEMP_DIR}/Run${2}.vcf.gz 2>>${log}
		tabix -p vcf ${TEMP_DIR}/Run${2}.vcf.gz 2>>${log}
	# I need to check if the tvc_out folder is still there.
	elif [ "`find ${1}/tvc*_out -maxdepth 0 -type d 2>/dev/null`" ]; then
		#	use the files already here.
		# the vcf files are not here. Exit with a file not found error. # Maybe I could call the run TVC script instead?
		checkFiles "${1}/tvc*_out/TSVC_variants.vcf" 
		cp ${1}/tvc*_out/TSVC_variants.vcf.gz ${TEMP_DIR}/Run${2}.vcf.gz
		cp ${1}/tvc*_out/TSVC_variants.vcf.gz.tbi ${TEMP_DIR}/Run${2}.vcf.gz.tbi
	else
		checkFiles "${1}/*.vcf"
		# what if there are multiple vcf files in the dir? Just take the first one it finds for now.
		vcf=`find ${1}/*.vcf -type f -printf "%f\n" | head -n 1`
		#Compress the vcf files and index them
		bgzip -c ${1}/${vcf} > ${TEMP_DIR}/Run${2}.vcf.gz 2>>${log}
		tabix -p vcf ${TEMP_DIR}/Run${2}.vcf.gz 2>>${log}
	fi
}

# gets only the chromosome and positions from the bed file to make bedtools coverage faster and the ouptut organized how we want it.
function setupBED {
	#python ${SOFTWARE_ROOT}/scripts/Remove_Overlap.py $BED 
	if [ "$CHR" != "" ]; then
		# get only the specified chromosome, and cut the other info so bedtools coverage will run faster , the output files will be smaller, and the number of columns in the output file will match for every project
		bed_name=`basename $1`
		grep -P "^${CHR}\t" $1 | bedtools sort -i | bedtools merge -i > ${TEMP_DIR}/${CHR}_sorted_merged_${bed_name}
		intersected_bed="${TEMP_DIR}/${CHR}_sorted_merged_${bed_name}"
	else
		# sort and merge the bed file
		bed_name=`basename $1`
		cat $1 | bedtools sort -i | bedtools merge -i > ${TEMP_DIR}/sorted_merged_${bed_name}
		intersected_bed="${TEMP_DIR}/sorted_merged_${bed_name}"
	fi
}

# $1: RUN1_DIR $2: run_num
function subset_low_cov {
	if [ "`find ${1}/cov_full -maxdepth 0 -type d 2>/dev/null`" ]; then
		checkFiles "${1}/cov_full/*.amplicon.cov.xls"
		amp="${1}/cov_full/*.amplicon.cov.xls"
	else
		checkFiles "${1}/*.amplicon.cov.xls"
		amp=`find ${1}/*.amplicon.cov.xls -type f`
	fi
	# Remove the header, and then pipe that into awk to get the amplicons with only > 30x coverage. Write it to a BED file.
	tail -n +2 $amp | \
		awk -v cutoff="$AMP_COV_CUTOFF" 'BEGIN{FS=OFS="\t";} {if ($10 >= cutoff) { printf("%s\t%s\t%s\t%s\t%s\n", $1, $2 - 1, $3, $4, $5); }}' \
		> ${TEMP_DIR}/run${2}_subset.bed
}

# $1: the bam file to index
function checkBamIndex {
	if [ ! "`find ${1}.bai -maxdepth 0 2>/dev/null`" ]; then
		java -jar $BAM_INDEXER INPUT=${1} OUTPUT=${1}.bai >>$log 2>&1 #2>&1 #hides the output of picard-tools.
	fi	
}

# Run GATK depth of coverage on a PTRIM.bam file in order to get total bases covered.
# It looks like Bedtools uses mapped and unmapped reads whereas Samtools and GATK use only mapped reads. Go to this site for more discovery: https://www.biostars.org/p/67579/
# $1: the PTRIM_RUN_BAM file, $2: the Overlap_subset.bed file $3: the run_number
function getDepths {
	bam_file=$1
	samtools depth \
		-b $2 \
		$bam_file \
		> ${TEMP_DIR}/Run$3_depths \
		&
#	bam_file=$1
#	checkBamIndex $bam_file
#	mkdir ${TEMP_DIR}/gatk${3}_out 2>/dev/null
#	java -jar $GATK \
#		--input_file $bam_file \
#		--analysis_type DepthOfCoverage \
#		--reference_sequence $REF_FASTA \
#		--intervals $2 \
#		-o ${TEMP_DIR}/gatk${3}_out/Run${3}_depths \
#		>> $log 2>&1 \
#		&
#	# bedtools does not need the .bam file to be indexed.
#	bedtools coverage \
#		-abam $bam_file \
#		-b $2 \
#		-d \
#		> ${TEMP_DIR}/Run${3}_depths \
#		&
}

# $1 is the bam file to use to run TVC, $2: the JSON PARAMETERS file used for TVC. $3: is the run number, $4: the location to output the PTRIM.bam
function runTVC {
	checkBamIndex $1
	mkdir ${TEMP_DIR}/tvc${3}_out 2>/dev/null
	# If the PTRIM.bam is already available, then don't generate it to save time.
	if [ "`find $4 -maxdepth 0 -type f 2>/dev/null`" ]; then
		# TVC uses the REF_FASTA file to call the variants first, and then intersects that with the --region-bed specified. Using the intersected_bed here would make no time saving difference
		${VARIANT_CALLER_DIR}/variant_caller_pipeline.py  \
			--input-bam $1 \
			--reference $REF_FASTA \
			--output-dir ${TEMP_DIR}/tvc${3}_out \
			--region-bed $BED \
			--parameters-file $2 \
			--hotspot-vcf ${TEMP_DIR}/final_Hotspot.vcf \
			--primer-trim-bed ${BED} \
			--bin-dir ${VARIANT_CALLER_DIR}  \
			>> $log 2>&1 \
			&
	# Generate the PTRIM.bam as it will be used to run samtools in this script
	else
		# TVC uses the REF_FASTA file to call the variants first, and then intersects that with the --region-bed specified. Using the intersected_bed here would make no time saving difference
		${VARIANT_CALLER_DIR}/variant_caller_pipeline.py  \
			--input-bam $1 \
			--reference $REF_FASTA \
			--output-dir ${TEMP_DIR}/tvc${3}_out \
			--region-bed $BED \
			--parameters-file $2 \
			--hotspot-vcf ${TEMP_DIR}/final_Hotspot.vcf \
			--postprocessed-bam=${4} \
			--primer-trim-bed ${BED} \
			--bin-dir ${VARIANT_CALLER_DIR}  \
			>> $log 2>&1 \
			&
	fi
}

# $1: The name of the program running
function waitForJobsToFinish {
	# Wait for job to finish. If there is an error, Fail will be 1 and the error message will be displayed 
	FAIL=0
	for job in `jobs -p`; do
		wait $job || let "FAIL+=1"
	done
	if [ "$FAIL" == "0" ]; then
		echo "	$OUTPUT_DIR $1 finished without problems at: `date`" 
		echo "$1 finished without problems at: `date`" >> $log
	else
		echo "-- ERROR: $1 had a problem. See ${log} for details --" 1>&2
		echo "-- ERROR: $1 had a problem. See ${log} for details --" >>$log
		exit 1
	fi
}


#for arguments
if [ $# -lt 15 ]; then # Technically there should be more than 20 arguments specified by the user, but more useful error messages can be displayed below
	usage
fi

RUNNING="Running QC_2Runs.sh. Option provided: "
counter=0
while :
do
	let "counter+=1"
	# If not enough inputs were given for an option, the while loop will just keep going. Stop it and print this error if it loops more than 100 times
	if [ $counter -gt 100 ]; then
		echo "USAGE: not all required inputs were given for options." 1>&2
		echo "$RUNNING"
		exit 8
	fi
	case $1 in
		-h | --help)
			usage
			;;
		-r | --run_dirs)
			RUN1_DIR=$2
			RUN2_DIR=$3
			RUNNING="$RUNNING --run_dirs: $2 $3 "
			shift 3
			;;
		-o | --output_dir)
			OUTPUT_DIR=$2
			RUNNING="$RUNNING --output_dir: $2 "
			shift 2
			;;
		-jo | --json_out)
			JSON_OUT=$2
			RUNNING="$RUNNING --json_out: $2 "
			shift 2
			;;
		-b | --project_bed)
			BED=$2
			RUNNING="$RUNNING --project_bed: $2 "
			shift 2
			;;
		-a | --amp_cov_cutoff)
			AMP_COV_CUTOFF=$2
			RUNNING="$RUNNING --amp_cov_cutoff: $2 "
			shift 2
			;;
		-jp | --json_paras)
			JSON_PARAS1=$2
			JSON_PARAS2=$3
			RUNNING="$RUNNING --json_paras: Json_paras1: $2, Json_paras2: $3 "
			shift 3
			;;
		-d | --depth_cutoffs) 
			DEPTH_CUTOFF1=$2
			DEPTH_CUTOFF2=$3
			RUNNING="$RUNNING --depth_cutoffs: Depth_cutoff1: $2, Depth_cutoff2: $3 "
			shift 3
			;;
		-gt | --gt_cutoffs)
			WT_CUTOFF1=$2
			HOM_CUTOFF1=$3
			WT_CUTOFF2=$4
			HOM_CUTOFF2=$5
			RUNNING="$RUNNING --gt_cutoffs: WT_1:$2 HOM_1:$3 WT_2:$4 HOM_2:$5 "
			shift 5
			;;

		-sb | --subset_bed)
			SUBSET_BED=$2
			RUNNING="$RUNNING --subset_bed: $2 "
			shift 2
			;;
		-v | --all_vcfs_dir)
			ALL_OUTPUT_DIR=$2
			RUNNING="$RUNNING --all_vcfs_dir: $2 "
			shift 2
			;;
		-chr | --subset_chr)
			CHR=$2
			RUNNING="$RUNNING --subset_chr: $2 "
			shift 2
			;;
		-B | --bases) 
			BASES="True"
			total_eligible_bases=$2
			total_possible_bases=$3
			RUNNING="$RUNNING --bases: total_eligible_bases: $2, total_possible_bases: $3 "
			shift 3
			;;
		-cl | --cleanup)
			CLEANUP="True"
			RUNNING="$RUNNING --cleanup "
			shift
			;;
		-*)
			printf >&2 'WARNING: Unknown option (ignored): %s\n' "$1"
			shift
			;;
		*)  # no more options. Stop while loop
			if [ "$1" != "" ]; then
				printf >&2 'WARNING: Unknown argument (ignored): %s\n' "$1"
				shift
			else
				break
			fi
			;;
	esac
done

# RUNNING contains all of the specified options
#echo "$RUNNING at `date`"

# Check to make sure the files actually exist
files=("$RUN1_DIR" "$RUN2_DIR" "${RUN1_DIR}/*.bam" "${RUN2_DIR}/*.bam" \
	"$BED" "$CDS_BED" "$SUBSET_BED" "$JSON_PARAS1" "$JSON_PARAS2" "$REF_FASTA" "$BAM_INDEXER" "$VARIANT_CALLER_DIR")
checkFiles $files

# ------------------------------------------------------------------------
# ------------ IF THERE ARE NO FILE ERRORS, START HERE -------------------
# ------------------------------------------------------------------------

mkdir -p $OUTPUT_DIR 2>/dev/null
TEMP_DIR="${OUTPUT_DIR}/temp_files"
mkdir $TEMP_DIR 2>/dev/null #This directory will hold all of the temporary files
log="${OUTPUT_DIR}/QC_2Runs.log"


# If the json files exist, use them. If not, make a new one
if [ "`find ${RUN1_DIR}/*.json* -maxdepth 0 -type f 2>/dev/null`" ]; then
	JSON1=`find ${RUN1_DIR}/*.json* -maxdepth 0 -type f | head -n 1`
else
	JSON1="${RUN1_DIR}/run_info.json_read"
fi
if [ "`find ${RUN2_DIR}/*.json* -maxdepth 0 -type f 2>/dev/null`" ]; then
	JSON2=`find ${RUN2_DIR}/*.json* -maxdepth 0 -type f | head -n 1`
else
	JSON2="${RUN2_DIR}/run_info.json_read"
fi


# If QC has already been run using the options specified, then no need to re-run it.
#qc_name=`basename $OUTPUT_DIR`
#if [ "`grep $qc_name $JSON_OUT 2>/dev/null`" ]; then
#	echo "	$OUTPUT_DIR has already been QCd. Skipping QC_2Runs.sh"
#else
# If the VCFs and depths have already been generated, then skip this step IF the variable $total_possible_bases and total_eligible_bases are passed in 
if [ "`find ${OUTPUT_DIR}/VCF1${CHR}_Final.vcf -type f 2>/dev/null`" -a \
	"`find ${OUTPUT_DIR}/VCF2${CHR}_Final.vcf -type f 2>/dev/null`" -a \
	"`find ${OUTPUT_DIR}/Both_Runs_${CHR}depths -maxdepth 0 -type f 2>/dev/null`" -a \
	"$BASES" == "True" ]; then
	echo "	$OUTPUT_DIR has already been QCd. Skipping QC_2Runs.sh"
else
	echo "	$OUTPUT_DIR Creating QC tables at `date`"
	#echo "For progress and results, see $log"
	echo "Creating QC tables for $OUTPUT_DIR at `date`" >>$log
	echo "----------------------------------------------" >>$log

	# this function finds the bam file for each run and subsets out the chr if specified 
	# The RUN1_BAM and RUN2_BAM variables are set up in this function
	setupBAMs $RUN1_DIR $RUN2_DIR

	# this function finds the vcf files and subsets out the chr if specified.
	setupVCF $RUN1_DIR 1
	setupVCF $RUN2_DIR 2
		
	# Get the unique and common variants to each vcf file.
	# I'm not sure if I actually have to do this
	echo "--- vcf-isec (column name warnings are fine) ---" >>${log}
	vcf-isec -f -c ${TEMP_DIR}/Run1.vcf.gz ${TEMP_DIR}/Run2.vcf.gz > ${TEMP_DIR}/VCF1_unique.vcf 2>>${log}
	vcf-isec -f -c ${TEMP_DIR}/Run2.vcf.gz ${TEMP_DIR}/Run1.vcf.gz > ${TEMP_DIR}/VCF2_unique.vcf 2>>${log}
	vcf-isec -f ${TEMP_DIR}/Run1.vcf.gz ${TEMP_DIR}/Run2.vcf.gz > ${TEMP_DIR}/VCF1_VCF2_common.vcf 2>>${log}

	# Merge the two compressed VCF files and remove the duplicates
	echo "--- vcf-merge ---" >>${log}
	vcf-merge ${TEMP_DIR}/Run1.vcf.gz ${TEMP_DIR}/Run2.vcf.gz --remove-duplicates > ${TEMP_DIR}/merged.vcf 2>>${log}

	# Now setup the BED files
	bed_name=`basename $BED`
	#bedtools intersect -a ${TEMP_DIR}/${intersected_bed} -b ${BED} -u -f 0.99 > ${TEMP_DIR}/intersect_${bed_name} 2>>$log
	intersected_bed="$BED"

	## interstect the subset.bed file with the specified subset bed files
	## if GET_CDS_DEPTHS is true, then the cds bed will not be subset out here so that all of the variants will be available, and  samtools depth can be run twice.
	#if [ "$CDS_BED" != "" -a "$GET_CDS_DEPTHS" != "True" ]; then
	#	cds_name=`basename $CDS_BED`
	#	# -f .99 option is not used here because the begin and end pos of the CDS region will not match up with the project bed file (it has intronic regions)
	#	bedtools intersect -a $BED -b ${CDS_BED} > ${TEMP_DIR}/intersect_${cds_name} 2>>$log
	#	intersected_bed="${TEMP_DIR}/intersect_${cds_name}"
	#fi
	if [ "$SUBSET_BED" != "" ]; then
		subset_name=`basename $SUBSET_BED`
		bedtools intersect -a ${intersected_bed} -b ${SUBSET_BED} -u -f 0.99 > ${TEMP_DIR}/intersect2_${subset_name} 2>>$log
		intersected_bed="${TEMP_DIR}/intersect2_${subset_name}"
	fi

	# This fuction gets only the chr and pos from the project bed file to calculate the total possible bases. (Also to save time and make the output better for bedtools coverage).
	setupBED ${intersected_bed} # The new bed file is also stored in intersected_bed
	# Calculate the total_possible_bases in the intersected bed file
	total_possible_bases=`awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $intersected_bed`

	# Function to keep only the amplicons that have depth coverage >= AMP COV CUTOFF.
	subset_low_cov $RUN1_DIR 1
	subset_low_cov $RUN2_DIR 2

	echo "--- Creating subset bed files (where amp cov > $AMP_COV_CUTOFF) ---" >>${log}
	# Now intersect the two subset bed files output by the subset_low_cov function
	bedtools intersect -a ${TEMP_DIR}/run1_subset.bed -b ${TEMP_DIR}/run2_subset.bed -u -f 0.99 > ${TEMP_DIR}/low_cov_subset.bed 2>>$log

	bed_name=`basename $intersected_bed`
	# Intersect the low_cov_subset.bed with the first bed specified. This shouldn't normally be needed, but it is a good precaution.
	bedtools intersect -a ${TEMP_DIR}/low_cov_subset.bed -b ${intersected_bed} -u -f 0.99 | bedtools sort -i | bedtools merge -i > ${TEMP_DIR}/subset_${bed_name} 2>>$log
	intersected_bed="${TEMP_DIR}/subset_${bed_name}"

	# check if we can subset the all comparison here
	if [ "$SUBSET_BED" != "" -a "$ALL_OUTPUT_DIR" != "" -a \
		"`find ${ALL_OUTPUT_DIR}/VCF1_Final.vcf -type f 2>/dev/null`" -a \
		"`find ${ALL_OUTPUT_DIR}/VCF2_Final.vcf -type f 2>/dev/null`" ]; then
		# subset out the SUBSET_BED file from the allVCF comparison
		echo "	Subsetting out the $intersected_bed file from the $ALL_OUTPUT_DIR comparison"
		echo "	Subsetting out the $intersected_bed file from the $ALL_OUTPUT_DIR comparison" >>$log
		bedtools intersect -header -a ${ALL_OUTPUT_DIR}/VCF1_Final.vcf -b $intersected_bed > ${OUTPUT_DIR}/VCF1_Final.vcf
		bedtools intersect -header -a ${ALL_OUTPUT_DIR}/VCF2_Final.vcf -b $intersected_bed > ${OUTPUT_DIR}/VCF2_Final.vcf
		# now subset the Depths files. Had to use python because the Both_Runs_depths is not in bed format.
		python ${QC_SCRIPTS}/subset_depths.py -b $SUBSET_BED -d ${ALL_OUTPUT_DIR}/Both_Runs_depths -o ${OUTPUT_DIR}/Both_Runs_depths
	fi

	# Only run this step if the Final VCF files have not already been created
	if [ ! "`find ${OUTPUT_DIR}/VCF1${CHR}_Final.vcf -type f 2>/dev/null`" -o ! "`find ${OUTPUT_DIR}/VCF2${CHR}_Final.vcf -type f 2>/dev/null`" ]; then
		# And then intersect that bed file with the merged vcf file to get only the variants that have > 30x coverage and that are found in the bed files specified.
		bedtools intersect -a ${TEMP_DIR}/merged.vcf -b ${intersected_bed} > ${TEMP_DIR}/merged_intersect.vcf 2>>$log
		
		# Filtering the Hotspot before vs after did not make a difference. Will filter before.
		# Remove the variants that have a multi-allelic call (i.e. A,G). Bonnie thinks they are a sequencing artifact.
		# if FAO + FRO is < Depth_Cutoff, that variant is removed
		echo "--- Filtering Variants (multi-allelic calls, and where FAO+FRO is < $DEPTH_CUTOFF1 in VCF1 and < $DEPTH_CUTOFF2 in VCF2 ---" >>${log}
		python ${QC_SCRIPTS}/QC_Filter.py \
			-v ${TEMP_DIR}/merged_intersect.vcf \
			-o ${TEMP_DIR}/filtered_merged_intersect.vcf \
			-O ${TEMP_DIR}/merged_variants_filtered_out.vcf \
			-c $DEPTH_CUTOFF1 -c $DEPTH_CUTOFF2 \
			>> $log 2>&1
		if [ "$?" != "0" ]
		then
			echo "ERROR: $OUTPUT_DIR had a problem filtering variants with QC_Filter.py... See $log for details" 1>&2
			exit 1
		fi
		
		echo "--- Creating the hotspot file ---" >>${log}
		# Create the hotspot file that has the variants for both of the runs.
		tvcutils prepare_hotspots -v ${TEMP_DIR}/filtered_merged_intersect.vcf -o ${TEMP_DIR}/final_Hotspot.vcf -r $REF_FASTA -s on -a on >>$log
		
		bgzip -c ${TEMP_DIR}/final_Hotspot.vcf > ${TEMP_DIR}/final_Hotspot.vcf.gz
		tabix -p vcf ${TEMP_DIR}/final_Hotspot.vcf.gz

		#The Hotspot file should be good to go now. Run TVC on the two original BAM files now using the awesome Hotspot file we created!
		echo "Running TVC using the generated hotspot file" >>$log
		echo "----------------------------------------------" >>$log
		# if CHR is nothing, then the PTRIM.bam will be written to the normal run's dir
		runTVC $RUN1_BAM $JSON_PARAS1 1${CHR} ${RUN1_DIR}/${CHR}/${CHR}PTRIM.bam
		runTVC $RUN2_BAM $JSON_PARAS2 2${CHR} ${RUN2_DIR}/${CHR}/${CHR}PTRIM.bam
		# after running the TVC result above, the TSVC_variants.vcf file created will include not only hotspot calls, but pretty much everything we filtered out prior to generating the final hotspot file. 
		#because TVC is actually run twice (1) as though there were no hotspot files defined (2) only using hotspot file. The final TSVC_variants.vcf reported is a combination of (1) and (2)  
		
		# wait for TVC to finish. Exit if TVC has a problem
		waitForJobsToFinish "TVC ${CHR}"

		# Now, intersect the variants to the hotspot file used as input in generating the final round of TSVC_variants.vcf. 
		# This way, we will match the no of variants in both runs and we will be able to proceed with generating the QC table (if the no of variants differed in the 2 runs, we would not be able to compare apples to apples)
		# get only the variants that are listed in the Hotspot file
		echo "--- vcf-isec to get only variants listed in the Hotspot file (column name warnings are fine) ---" >>${log}
		vcf-isec -f ${TEMP_DIR}/tvc1${CHR}_out/TSVC_variants.vcf.gz ${TEMP_DIR}/final_Hotspot.vcf.gz > ${TEMP_DIR}/VCF1_Intersect_Hotspot.vcf 2>>${log}
		vcf-isec -f ${TEMP_DIR}/tvc2${CHR}_out/TSVC_variants.vcf.gz ${TEMP_DIR}/final_Hotspot.vcf.gz > ${TEMP_DIR}/VCF2_Intersect_Hotspot.vcf 2>>${log}
		
		# Finally, if the user specified a subset chr, put that in the name of the final VCF file.
		VCF1_FINAL="${OUTPUT_DIR}/VCF1${CHR}_Final.vcf"
		VCF2_FINAL="${OUTPUT_DIR}/VCF2${CHR}_Final.vcf"
		echo "--- removing duplicates entries---" >>${log}
		# I keep trying to think of ways to shortcut around having to do this extra filtering, but Ozlem is doing it all
		# because there were special cases where these steps were necessary. So I'll keep them.
		# We need to remove the duplicate entries again.
		python ${QC_SCRIPTS}/QC_Filter.py \
			-v ${TEMP_DIR}/VCF1_Intersect_Hotspot.vcf \
			-o $VCF1_FINAL \
			-c $DEPTH_CUTOFF1 \
			>> $log 2>&1
		python ${QC_SCRIPTS}/QC_Filter.py \
			-v ${TEMP_DIR}/VCF2_Intersect_Hotspot.vcf \
			-o $VCF2_FINAL \
			-c $DEPTH_CUTOFF2 \
			>> $log 2>&1
	fi

	# If the depths have already been generated, then skip this step
	if [ ! "`find ${OUTPUT_DIR}/Both_Runs_${CHR}depths -maxdepth 0 -type f 2>/dev/null`" ]; then
		# get Depths. Not needed for hotspot creation and such, but it will be used to find the total # of eligible bases later on.
		echo "${RUN1_DIR}/${CHR}/${CHR}PTRIM.bam" > ${TEMP_DIR}/depths_files 
		echo "${RUN2_DIR}/${CHR}/${CHR}PTRIM.bam" >> ${TEMP_DIR}/depths_files 

		# Use samtools depth to get the depths of both Run's ptrim.bam files at every position listed in the intersected_bed.
		# -f option requires a file with a bam file on each line in that file
		samtools depth \
			-f ${TEMP_DIR}/depths_files \
			-b $intersected_bed \
			> ${TEMP_DIR}/Both_Runs_${CHR}depths \
			&

		waitForJobsToFinish "samtools depths ${CHR}"
		# Move the depths to the Output dir in case the script is interrupted in the middle of running samtools
		mv  ${TEMP_DIR}/Both_Runs_${CHR}depths ${OUTPUT_DIR}

	fi
	# Use awk to get the total_eligible_bases by getting only the base positions from both samtools depth outputs that have greater than the cutoff depth.
	total_eligible_bases=`awk -v cutoff1=$DEPTH_CUTOFF1 -v cutoff2=$DEPTH_CUTOFF2 '{ if ($3 >= cutoff1 && $4 >= cutoff2) printf "."}' ${OUTPUT_DIR}/Both_Runs_${CHR}depths | wc -c`
fi

# This script takes the two VCF files generated from running TVC using the same hotspot file, matches them and outputs the info to the csv and the json_out file.
 # total_possible_bases is calculated after the  setupBed function.
  # The ouptut json file will be used to generate the 3x3 QC tables.
python ${QC_SCRIPTS}/QC_Compare_VCFs.py \
	--vcfs ${OUTPUT_DIR}/VCF1${CHR}_Final.vcf ${OUTPUT_DIR}/VCF2${CHR}_Final.vcf \
	--jsons $JSON1 $JSON2 \
	--gt_cutoffs $WT_CUTOFF1 $HOM_CUTOFF1 $WT_CUTOFF2 $HOM_CUTOFF2 \
	--bases $total_eligible_bases $total_possible_bases \
	--out_csv ${OUTPUT_DIR}/matched_variants${CHR}.csv \
	--json_out ${JSON_OUT} \
	--chr "$CHR" # This option will add a prefix of the chromosome used. If $CHR is blank or nothing, no prefix will be added. That way both chr1 and all exome data can be run together without any problems That way both chr1 and all exome data can be run together without any problems
if [ $? -ne 0 ]; then
	echo "ERROR: $OUTPUT_DIR QC_Compare_VCFs.py had a problem!! " 1>&2
	exit 
fi

# Cleanup is done in cleanup.py called by QC_sample
# Cleanup and done.
#if [ "$CLEANUP" == "True" ]; then
#	rm -rf ${TEMP_DIR}
#	# Remove the chr subsets if they were created.
#	if [ "$CHR" != "" ]; then
#		rm -rf ${RUN1_DIR}/${CHR} ${RUN2_DIR}/${CHR} 2>/dev/null
#	fi
#	# Remove the PTRIM.bam
#	rm ${RUN1_DIR}/PTRIM.bam* ${RUN2_DIR}/PTRIM.bam* 2>/dev/null
#fi

echo "	$OUTPUT_DIR Finished QC." 
#echo "----------------------------------------------"
echo "$OUTPUT_DIR Finished QC." >>$log
#fi

