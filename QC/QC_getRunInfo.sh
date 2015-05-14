#!/bin/bash

# Goal: Run GATK on the given sample using the three different BED files to determine the quality of coverage for this sample.
# Also gather and generate the other Run info from the backupPDF.pdf file or the TS API
# TS version: 4.2
# Author: Jeff L

# QSUB variables

#$ -S /bin/bash
#$ -cwd
#$ -N QC
#$ -o qc_sample_out.log
#$ -e qc_sample_out.log
#$ -q plugin.q ##TEMP FOR WALES
#$ -V


# Define file paths
SOFTWARE_ROOT="/rawdata/legos"
QC_SCRIPTS="${SOFTWARE_ROOT}/scripts/QC"
BAM_INDEXER='/opt/picard/picard-tools-current/BuildBamIndex.jar'
GATK="/results/plugins/variantCaller/TVC/jar/GenomeAnalysisTK.jar"
REF_FASTA="/results/referenceLibrary/tmap-f3/hg19/hg19.fasta"
VARIANT_CALLER_DIR="/results/plugins/variantCaller"
TVC_VERSION='4.2'

function usage {
cat << EOF
USAGE: bash QC_getRunInfo.sh 
	-h | --help
	-r | --run_dir <path/to/run_dir> (run dir should have the vcf file, the amplicon.cov.xls file, and the .json file)
	-o | --out_dir <path/to/output_dir> 
	-a | --amp_cov_cutoff <min_amp_coverage> (Cutoff for # of amplicon reads.) 
	-d | --depth_cutoff <min_base_depth> (Cutoff for base depth)
	-wh | --wt_hom_cutoff <WT_Cutoff> <HOM_Cutoff> (WT and HOM cutoffs)
	-pb | --project_bed <path/to/project.bed	(To generate the PTRIM.bam. TEMP FOR WALES)
	-pj | --ptrim_json <path/to/tvc_json>	(Regenerate the PTRIM.bam if it is not found. Used to calculate depths.)
	-bb | --beg_bed <path/to/beginning_loci_bed> (Only available if the Run has a PTRIM.bam. Run GATK using the bed file with only the 10th pos of the amplicons)
	-eb | --end_bed <path/to/end_loci_bed> (Only available if the Run has a PTRIM.bam. Run GATK using the bed file with only the 10th pos from the end of the amplicons)
	-cb | --cds_bed <path/to/CDS_bed> (Optional. Only available if the Run has a PTRIM.bam. Run GATK on the CDS region of the bed file)
	-pl | --pool_dropout	(Optional. Will include the pool dropout script.)
	-cl | --cleanup (Optional. Will delete the file specified as --out_dir after generating the QC metrics needed.)
EOF
exit 8
}


# Checks to ensure that the files provided exist and are readable. 
# @param $1 should be a list of files.
function checkFiles {
	files=$1
	for file in ${files[@]}; do
		if [ ! "`find ${file} -maxdepth 0 2>/dev/null`" ]; then
			echo "-- ERROR -- '${file}' not found or not readable" 1>&2
			exit 4
		fi
	done
}

# $1: The name of the program running
function waitForJobsToFinish {
	# Wait for job to finish. If there is an error, Fail will be 1 and the error message will be displayed 
	FAIL=0
	for job in `jobs -p`; do
		wait $job || let "FAIL+=1"
	done
	if [ "$FAIL" == "0" ]; then
		if [ "$job" != "" ]; then
			echo "	$1 finished successfully at: `date`"
		fi
	else
		echo "--- ERROR: $1 had a problem. Not copying data ---" 1>&2
		# Don't exit because we can still get the variant info.
		#exit 1
	fi
}

# For arguments
if [ $# -lt 5 ]; then
	usage
fi

RUNNING="Running QC_getRunInfo.sh. Options provided: "
while :
do
	case $1 in
		-h | --help)
			usage
			;;
		-r | --run_dir)
			RUN_DIR=$2
			RUNNING="$RUNNING --run_dir: $2 "
			shift 2
			;;
		-o | --out_dir)
			OUTPUT_DIR=$2
			RUNNING="$RUNNING --out_dir: $2 "
			shift 2
			;;
		-a | --amp_cov_cutoff)
			AMP_COV_CUTOFF=$2
			RUNNING="$RUNNING --amp_cov_cutoff: $2 "
			shift 2
			;;
		-d | --depth_cutoff)
			DEPTH_CUTOFF=$2
			RUNNING="$RUNNING --depth_cutoff: $2 "
			shift 2
			;;
		-wh | --wt_hom_cutoff)
			WT_CUTOFF=$2
			HOM_CUTOFF=$3
			RUNNING="$RUNNING --wt_hom_cutoff: WT: $2, HOM: $3 "
			shift 3
			;;
		-pb | --project_bed)
			PROJECT_BED=$2
			RUNNING="$RUNNING --project_bed: $2 "
			shift 2
			;;
		-pj | --ptrim_json)
			REGEN_PTRIM="True"
			TVC_JSON=$2
			RUNNING="$RUNNING --ptrim_json: $2 "
			shift 2
			;;
		-bb | --beg_bed)
			BEG_BED=$2
			RUNNING="$RUNNING --beg_bed: $2 "
			shift 2
			;;
		-eb | --end_bed)
			END_BED=$2
			RUNNING="$RUNNING --end_bed: $2 "
			shift 2
			;;
		-cb | --cds_bed)
			CDS_BED=$2
			RUNNING="$RUNNING --cds_bed: $2 "
			shift 2
			;;
		-pl | --pool_dropout)
			POOL_DROPOUT="True"
			RUNNING="$RUNNING --pool_dropout "
			shift
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
echo "$RUNNING at `date`"

# Check to ensure that all of the required options are provided
if [ "$RUN_DIR" == "" -o "$OUTPUT_DIR" == "" -o "$AMP_COV_CUTOFF" == "" -o "$DEPTH_CUTOFF" == "" \
	-o "$WT_CUTOFF" == "" -o "$HOM_CUTOFF" == "" ]; then
	echo "--run_dir, --output_dir, --amp_cov_cutoff, --depth_cutoff, --wt_hom_cutoff are all required options"
	exit 8
fi
# -o "$BEG_BED" == "" -o "$END_BED" == "" 
# , --begbed, --endbed 


files=("$RUN_DIR" "$BEG_BED" "$END_BED" "$REF_FASTA" "$GATK")
# Don't check for the .json file yet because if it doesn't exist, we can create it.
#files=("$RUN_DIR" ${RUN_DIR}/*.json* "$BEG_BED" "$END_BED" "$REF_FASTA" "$GATK")
checkFiles $files

PTRIM_BAM=`find ${RUN_DIR}/PTRIM.bam -maxdepth 0 -type f 2>/dev/null`

# get the VCF file
if [ "`find ${RUN_DIR}/tvc*_out -maxdepth 0 -type d 2>/dev/null`" ]; then
	checkFiles "${RUN_DIR}/tvc*_out/TSVC_variants.vcf" 
	VCF=`find ${RUN_DIR}/tvc*_out/TSVC_variants.vcf -type f`
else
	checkFiles "${RUN_DIR}/*.vcf"
	VCF=`find ${RUN_DIR}/*.vcf -type f | head -n 1`
fi
if [ "`find ${RUN_DIR}/cov_full -maxdepth 0 -type d 2>/dev/null`" ]; then
	checkFiles "${RUN_DIR}/cov_full/*.amplicon.cov.xls"
	AMP=`find ${RUN_DIR}/cov_full/*.amplicon.cov.xls -type f | head -n 1`
else
	checkFiles "${RUN_DIR}/*.amplicon.cov.xls"
	AMP=`find ${RUN_DIR}/*.amplicon.cov.xls -type f | head -n 1`
fi


mkdir -p $OUTPUT_DIR
log="${OUTPUT_DIR}/getRunInfo.log"
echo "$RUNNING at `date`" >> $log

# Check to see if the depths already exist
if [ "`find ${OUTPUT_DIR}/forward_beg_depths -type f 2>/dev/null`" -a "`find ${OUTPUT_DIR}/reverse_beg_depths -type f 2>/dev/null`" -a \
		"`find ${OUTPUT_DIR}/forward_end_depths -type f 2>/dev/null`" -a "`find ${OUTPUT_DIR}/reverse_end_depths -type f 2>/dev/null`" -a \
		"`find ${OUTPUT_DIR}/overall_run_depths -type f 2>/dev/null`" ]; then
	echo "	$OUTPUT_DIR already has the depth files. Skipping samtools depth and getting the metrics"
else
	# Check to see if the run has a PTRIM.bam. if it does, then we can run samtools depth.
	if [ "$PTRIM_BAM" == "" -a "$REGEN_PTRIM" == "True" ]; then
		echo "	$RUN_DIR has no PTRIM. Regenerating it now at: `date`"
	    bam=`find ${RUN_DIR}/*.bam -maxdepth 0 -type f 2>/dev/null | head -n 1`
	
		#now run TVC v4.2. 
		mkdir -p ${OUTPUT_DIR}/tvc${TVC_VERSION}_out 2>/dev/null
		${VARIANT_CALLER_DIR}/variant_caller_pipeline.py \
			--input-bam "$bam" \
			--reference-fasta "$REF_FASTA" \
			--output-dir "${OUTPUT_DIR}/tvc${TVC_VERSION}_out" \
			--parameters-file "$TVC_JSON" \
			--bin-dir "${VARIANT_CALLER_DIR}" \
			--region-bed  "${PROJECT_BED}" \
			--primer-trim-bed "${PROJECT_BED}" \
		    --postprocessed-bam=${RUN_DIR}/PTRIM.bam \
	    	> ${OUTPUT_DIR}/tvc${TVC_VERSION}_out/log.out 2>&1
	
		# check to make sure TVC finished correctly.
		if [ $? -ne 0 ]; then
			echo "	$RUN_DIR was unable to generate the PTRIM.bam at `date`"
			echo "	$RUN_DIR was unable to generate the PTRIM.bam at `date`" >> $log
		else
			PTRIM_BAM="${RUN_DIR}/PTRIM.bam"
			echo "	$PTRIM_BAM was generated successfully at `date`"
			echo "	$PTRIM_BAM was generated successfully at `date`" >> $log
		fi
	fi
	
	# 
	if [ "$PTRIM_BAM" == "" ]; then
		echo "	PTRIM.bam was not available for this run. Skipping samtools depth"
		echo "	PTRIM.bam was not available for this run. Skipping samtools depth" >>$log
	else
		echo "	$PTRIM_BAM beginning samtools depth at `date`"
		mkdir ${OUTPUT_DIR}/coverages 2>/dev/null
		# get the depth of the beginning and end of the forward and reversestrands
		# -g means include, and -G means exclude. 0x10 means reverse strands
		samtools depth -b $BEG_BED -G 0x10 $PTRIM_BAM > ${OUTPUT_DIR}/coverages/forward_beg_depths &
		samtools depth -b $BEG_BED -g 0x10 $PTRIM_BAM > ${OUTPUT_DIR}/coverages/reverse_beg_depths &
		samtools depth -b $END_BED -G 0x10 $PTRIM_BAM > ${OUTPUT_DIR}/coverages/forward_end_depths &
		samtools depth -b $END_BED -g 0x10 $PTRIM_BAM > ${OUTPUT_DIR}/coverages/reverse_end_depths &
		samtools depth -b $PROJECT_BED $PTRIM_BAM > ${OUTPUT_DIR}/coverages/overall_run_depths &
	fi
	
	waitForJobsToFinish "$PTRIM_BAM samtools depth"
	# Move the depths files to the output dir. This way, they will only be in the output dir if samtools depth finished successfully on everything
	mv ${OUTPUT_DIR}/coverages/* $OUTPUT_DIR 2>/dev/null
	rmdir ${OUTPUT_DIR}/coverages 2>/dev/null
fi	

# The depth files should already be generated by this point, so we can just get the stats we need here
# get the number of amplicons that have the 10th forward strand base covered at $AMP_COV_CUTOFF
forward_begbpCov=`awk -v cutoff=$AMP_COV_CUTOFF '{ if($3 >= cutoff) printf "."}' ${OUTPUT_DIR}/forward_beg_depths | wc -c 2>/dev/null`
# get the number of amplicons that have the 10th reverse strand base covered at $AMP_COV_CUTOFF
reverse_begbpCov=`awk -v cutoff=$AMP_COV_CUTOFF '{ if($3 >= cutoff) printf "."}' ${OUTPUT_DIR}/reverse_beg_depths | wc -c 2>/dev/null`
# get the number of amplicons that have the 10th forward strand base from the end covered at $AMP_COV_CUTOFF
forward_endbpCov=`awk -v cutoff=$AMP_COV_CUTOFF '{ if($3 >= cutoff) printf "."}' ${OUTPUT_DIR}/forward_end_depths | wc -c 2>/dev/null`
# get the number of amplicons that have the 10th reverse strand base from the end covered at $AMP_COV_CUTOFF
reverse_endbpCov=`awk -v cutoff=$AMP_COV_CUTOFF '{ if($3 >= cutoff) printf "."}' ${OUTPUT_DIR}/reverse_end_depths | wc -c 2>/dev/null`
# get the number of bases covered at $DEPTH_CUTOFF
base_cov=`awk -v cutoff=$DEPTH_CUTOFF '{ if($3 >= cutoff) printf "."}' ${OUTPUT_DIR}/overall_run_depths | wc -c 2>/dev/null`

# Calculate the total_possible_bases in the project bed file. NOT the project bed file, the merged bed file!
total_possible_bases=`cat $PROJECT_BED | bedtools sort -i | bedtools merge -i | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'`

num_amps=`tail -n +2 $AMP | wc -l`

# add up the beginning of the forward reads depths with the end of the reverse reads depths to get the "30x coverage at beginning of amplicon"
begin_amp_cov=$(( (forward_begbpCov + reverse_endbpCov) / 2 ))
#echo "$begin_amp_cov"
# add up the end of the forward reads depths with the beginning of the reverse reads depths to get the "30x coverage at end of amplicon"
end_amp_cov=$(( (forward_endbpCov + reverse_begbpCov) / 2 ))
#echo "$end_amp_cov"

# Remove the header, and then pipe that into awk to get the amplicons with only > 30x coverage.
amp_cov=`tail -n +2 $AMP | awk -v cutoff=$AMP_COV_CUTOFF '{ if ($10 >= cutoff) printf "."}' | wc -c` 

# Find this run's .json file, or add another one
if [ "`find ${RUN_DIR}/*.json* -maxdepth 0 2>/dev/null`" ]; then
	JSON=`find ${RUN_DIR}/*.json* -maxdepth 0`
else
	run_name=`basename $RUN_DIR`
	JSON="${RUN_DIR}/${run_name}.json_read"
fi

# filter the VCF file to then get the TS_TV ratio
python ${QC_SCRIPTS}/QC_Filter.py \
	-v $VCF \
	-o ${OUTPUT_DIR}/filtered${DEPTH_CUTOFF}x.vcf \
	-O ${OUTPUT_DIR}/leftover_variants${DEPTH_CUTOFF}x.vcf \
	-c $DEPTH_CUTOFF \
	-j $JSON >> $log
TS_TV=`cat ${OUTPUT_DIR}/filtered${DEPTH_CUTOFF}x.vcf | vcf-tstv | grep -Po "\d+\.\d+"`

# get the medianReadCoverageOverall for each amplicon by adding the total forward and reverse reads columns in the .amplicon.cov.xls, and then finding the median.
medainReadCoverageOverall=`awk '{ print ($11 + $12) }' $AMP | awk '{ count[NR] = $1; } END { if (NR % 2) { print count[(NR + 1) / 2]; } else { print (count[(NR / 2)] + count[(NR / 2) + 1]) /2.0; }}'`

# These QC metrics will be added to the json file of the run.
metrics="base_cov:total_possible_bases:num_amps:amp_cov:begin_amp_cov:end_amp_cov:ts_tv:median_coverage_overall;$base_cov:$total_possible_bases:$num_amps:${amp_cov}:$begin_amp_cov:$end_amp_cov:$TS_TV:$medainReadCoverageOverall" 

# Check if this run has a reportPDF.pdf
if [ "`find ${RUN_DIR}/report.pdf -type f 2>/dev/null`" ]; then
	# add the pdf file to QC_printRunInfo.py
	#python ${QC_SCRIPTS}/QC_parse_TS_PDF.py -j $JSON -p ${RUN_DIR}/backupPDF.pdf -o $OUTPUT_DIR >> $log 2>&1
	python ${QC_SCRIPTS}/QC_parse_TS_PDF.py -j $JSON -p ${RUN_DIR}/report.pdf -o $OUTPUT_DIR >> $log 2>&1
	if [ $? -ne 0 ]; then
		echo "	ERROR: QC_parse_TS_PDF.py was unsuccessful. See $OUTPUT_DIR/getRunInfo.log for details" 1>&2
		ERRORS="True"
	else
		echo "$RUN_DIR Gathered the pdf info at `date`" >>$log
	fi
# If no reportPDF is found, then check if there is a backupPDF. If not, hopefully these metrics were already gathered. 
elif [ "`find ${RUN_DIR}/backupPDF.pdf -type f 2>/dev/null`" ]; then
	python ${QC_SCRIPTS}/QC_parse_TS_PDF.py -j $JSON -p ${RUN_DIR}/backupPDF.pdf -o $OUTPUT_DIR >> $log 2>&1
    if [ $? -ne 0 ]; then
		echo "  ERROR: QC_parse_TS_PDF.py was unsuccessful. See $OUTPUT_DIR/getRunInfo.log for details" 1>&2
		ERRORS="True"
	else
		echo "$RUN_DIR Gathered the pdf info at `date`" >>$log
	fi
fi

#echo "$VCF is now getting var info"
# get the total number of variants in the vcf file, and print the other GATK metrics.
python ${QC_SCRIPTS}/QC_getVarInfo.py \
	--filtered_vcf ${OUTPUT_DIR}/filtered${DEPTH_CUTOFF}x.vcf \
	--json $JSON \
	--metrics "$metrics" \
	--gt_cutoffs $WT_CUTOFF $HOM_CUTOFF \
	>>$log 2>&1

if [ $? -ne 0 ]; then
	echo "	ERROR QC_getVarInfo.py was unsuccessful. See $log for details" 1>&2
	ERRORS="True"
#else
#	echo "$RUN_DIR Gathered the vcf info at `date`">> $log
fi

# parameter for pool dropout to only run on exon or otherwise specified data.
if [ "$POOL_DROPOUT" == "True" ]; then
	#run Matt's script to add pool info to the json file 
	echo "	Generating median read length info for each pool"
	python2.7  ${QC_SCRIPTS}/pool_dropout.py \
		-c ${OUTPUT_DIR}/reverse_beg_depths \
		-c ${OUTPUT_DIR}/forward_beg_depths \
		-j $JSON \
		-o ${OUTPUT_DIR} \
		-b "${PROJECT_BED}" \
		-p 10 
fi

# If there were errors, then exit with a status of 1
if [ "$ERRORS" == "True" ]; then
	echo "	$RUN_DIR finished with errors at `date`"
	echo "$RUN_DIR finished with errors at `date`" >>$log
	exit 1
fi

echo "	$RUN_DIR finished getting run info at `date`"
echo "$RUN_DIR finished getting run info at `date`" >>$log

# cleanup and finished
if [ "$CLEANUP" == "True" ]; then
	# Keep the depth files so that we can use them for later calculations
	rm -rf $OUTPUT_DIR/tvc${TVC_VERSION}_out
fi

