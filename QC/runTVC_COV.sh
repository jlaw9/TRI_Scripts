#! /bin/bash

# Goal: Run coverage analysis and TVC on a BAM file 
# It would be great to allow you to either specify a sample, plate, or project, or a csv.

# Parameters for Qsub:
#$ -S /bin/bash
#$ -N TVC_Cov
#$ -o tvc_results.log
#$ -e tvc_results.log
#$ -q all.q
#$ -cwd
#$ -V

# Default directories
SOFTWARE_ROOT="/rawdata/legos" # used for flags.sh
PICARD_TOOLS_DIR="/opt/picard/picard-tools-current"
COV_ANALYSIS_DIR="/results/plugins/coverageAnalysis"
VARIANT_CALLER_DIR="/results/plugins/variantCaller"
REF_FASTA="/results/referenceLibrary/tmap-f3/hg19/hg19.fasta"

# Default options
RUN_COV="False"
AMPLISEQ="True"
RUN_TVC="False"
TVC_VERSION='4.2'
CLEANUP="False"
NOERRS="True"
FORCED="False"

function usage {
cat << EOF
USAGE: bash runTVC_Cov.sh <path/to/BAM_File.bam> [OPTIONS]
Options:
	-A | --ampliseq 	(Default. Will run coverage analysis and TVC with Ampliseq options)
	-T | --targetseq	(Will run coverage analysis with Targetseq options)
	-c | --cov <path/to/Merged_BED> <path/to/Unmerged_BED>   (If no .cov.xls file is present in the bam file's directory, use this option to run Coverage Analysis)
	-t | --tvc <path/to/Project_BED> <path/to/tvc_json> 	(region bed file and parameter settings file which will be used by TVC.)
	-p | --ptrim <PTRIM.bam>	(Will generate the PTRIM.bam if it is not already found)
	-th | --tvc_hotspot <path/to/Hotspot.vcf> 	(The hotspot file used to run TVC)
	-o | --output_dir <path/to/Output_Dir>	(Output files will be placed here. Default is the Bam File's directory)
	-rdf | --remove_dup_flags	(Will use samtools to check if there are duplicate flags. If there are, remove the duplicate flags from the bam file using Picard tools)
	-fd | --flag_dups		(Will place the duplicate flag in the bam file)
	-cl | --cleanup		(Delete the temporary files generated to run TVC and Coverage Analysis. If there are any errors, the temp files are not deleted.)
	-f | --forced			(Use this option to run tvc/cov even if there is a .vcf/.cov.xls file already present)
EOF
exit 8
}

#for arguments
if [ $# -lt 3 ]; then
	usage
fi

RUNNING="Starting runTVC_Cov.sh with these options: "
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
		*.bam)
			BAM_FILE=$1
			RUNNING="$RUNNING BAM: $1 "
			shift 
			;;
		-A | --ampliseq)
			AMPLISEQ="True"
			RUNNING="$RUNNING --ampliseq "
			shift 
			;;
		-T | --targetseq)
			TARGETSEQ="True"
			AMPLISEQ="False"
			RUNNING="$RUNNING --targetseq "
			shift 
			;;
		-c | --cov)
			RUN_COV="True"
			MERGED_BED=$2
			UNMERGED_BED=$3
			RUNNING="$RUNNING --cov merged_bed: $2, unmerged_bed: $3 "
			shift 3
			;;
		-t | --tvc)
			RUN_TVC="True"
			TVC_BED=$2
			TVC_JSON=$3
			RUNNING="$RUNNING --tvc bed: $2 json: $3 "
			shift 3
			;;
		-p | --ptrim)
			PTRIM="$2"
			RUNNING="$RUNNING --ptrim: $2 "
			shift 2
			;;
		-th | --tvc_hotspot)
			TVC_HOTSPOT=$2
			RUNNING="$RUNNING --tvc_hotspot: $2 "
			shift 2
			;;
		-o | --output_dir)
			OUTPUT_DIR=$2
			RUNNING="$RUNNING --output_dir: $2 "
			shift 2
			;;

		-rdf | --remove_dup_flags)
			REMOVE_DUP_FLAGS="True"
			RUNNING="$RUNNING --remove_dup_flags "
			shift
			;;
		-fd | --flag_dups)
			FLAG_DUPS="True"
			RUNNING="$RUNNING --flag_dups "
			shift
			;;
		-cl | --cleanup)
			CLEANUP="True"
			RUNNING="$RUNNING --cleanup "
			shift
			;;
		-f | --forced)
			FORCED="True"
			RUNNING="$RUNNING --forced"
			shift
			;;
		-*)
			printf >&2 'WARNING: Unknown option (ignored): %s\n' "$1"
			shift
			;;
		 *)  # no more options. Stop while loop
			if [ "$1" != '' ]; then
					printf >&2 'WARNING: Unknown option (ignored): %s\n' "$1"
					shift
			else    
				break
			fi  
			;;  
	esac
done

# RUNNING is a variable containing all of the options specified
echo "$RUNNING at `date`"

# Now check the files to see if the bam file exists.
if [ "$BAM_FILE" == '' ]; then
	echo "ERROR: MUST SPECIFY A BAM FILE" 1>&2
#	usage
	exit 4
fi
if [ ! "`find $BAM_FILE -maxdepth 0 2>/dev/null`" ]; then
	echo "ERROR -- '${BAM_FILE}' not found or not readable" 1>&2
	exit 8
fi
# If no output dir is specified, default will be the bam file's directory
if [ "$OUTPUT_DIR" == '' ]; then
	OUTPUT_DIR=`dirname $BAM_FILE`
fi

# ------------------------------------------------
# ---------------- FUNCTIONS ---------------------
# ------------------------------------------------

function checkFlags {
	# Remove or flag the duplicates according to what the user specified.
	if [ "$REMOVE_DUP_FLAGS" == "True" ]; then
		bash ${SOFTWARE_ROOT}/scripts/flags.sh \
			$BAM_FILE \
			--remove_dup \
			--cleanup
	elif [ "$FLAG_DUPS" == "True" ]; then
		bash ${SOFTWARE_ROOT}/scripts/flags.sh \
			$BAM_FILE \
			--flag_dups 

		bam_dir=`dirname $BAM_FILE`
		bam_name=`basename $BAM_FILE`
		if [ "`find ${bam_dir}/Flag_${bam_name} -maxdepth 0 -type f 2>/dev/null`" ]; then
			mkdir ${bam_dir}/Bam_Backup
			mv ${BAM_FILE}* ${bam_dir}/Bam_Backup
			BAM_FILE="${bam_dir}/Flag_${bam_name}"
		fi
		#	--cleanup
	fi
}

# $1: the bam file to index, $2 the log file
function checkBamIndex {
	# If the bam file needs to be indexed, index the bam file
	if [ ! "`find ${1}.bai -maxdepth 0 2>/dev/null`" ]; then
		echo "	Indexing $BAM_FILE at `date`"
		java -jar ${PICARD_TOOLS_DIR}/BuildBamIndex.jar INPUT=${1} OUTPUT=${1}.bai >>$2 2>&1 

		# If for some reason the bam file couldn't be indexed, then quit. The indexed bam file is needed to Cov and TVC
		if [ $? -ne 0 ]; then
			echo "	First index failed... sorting and retrying with samtools: $BAM_FILE at `date`"
			bam_dir=`dirname $1`
			samtools sort $1 ${bam_dir}/sorted
			mv ${bam_dir}/sorted.bam $1
			samtools index $1
			if [ $? -ne 0 ]; then
				echo "ERROR: Unable to index $1 even after sorting with samtools sort ... quitting"
				exit 1
			fi
		fi
	fi	
}

# Runs coverage analysis, then copies the .amplicon.cov.xls file it needs, and deletes the other files generated.
function runCov {
	mkdir -p ${OUTPUT_DIR}/cov_full 2>/dev/null
	checkBamIndex $BAM_FILE ${OUTPUT_DIR}/cov_full/log.out

	echo "	$BAM_FILE beginning Coverage Analysis at: `date`"
	
	# -a is for Ampliseq.  -g option gives the gene name and GC content  -D specifies the output directory, -d is for bam files with flagged duplicates
	# Ozlem had to modify the run_coverage_analysis.sh to include the 30x coverage. She modified the targetReadStats.pl which is being called by run_cv_analysis.sh
	# For some reason, putting the merged BED file for -A and the unmerged BED file for -B with the -g option generated the same .amplicon.cov.xls file as the web browser did. We will use these settings.
	run_cov="${COV_ANALYSIS_DIR}/run_coverage_analysis.sh"
	if [ "$AMPLISEQ" == "True" ]; then
		run_cov="$run_cov -ag"
	elif [ "$TARGETSEQ" == "True" ]; then
		run_cov="$run_cov -c -d"
	fi
	run_cov="""$run_cov -D "${OUTPUT_DIR}/cov_full" -A "${MERGED_BED}" -B "${UNMERGED_BED}" "$REF_FASTA" "$BAM_FILE"""" 	

	echo "	running coverage analysis with these options: $run_cov at `date`" >> ${OUTPUT_DIR}/cov_full/log.out
	$run_cov >> ${OUTPUT_DIR}/cov_full/log.out 2>&1

	#wait for job to finish. If there is an error in coverage analysis, $? will be 1 and the error message will be displayed 
	if [ $? -eq 0 ]; then
		# Coverage analysis was successful. Copy the .amplicon.cov.xls file, and Delete the other temporaray files.
		echo "	$BAM_FILE Coverage analysis finished successfully at: `date`"
	else
		NOERRS="False"
		# Something went wrong with coverage analysis. Not copying the data.
		echo 1>&2
		echo "--- ERROR: $BAM_FILE Coverage analysis was unsuccessful. Not copying or deleting data ---" 1>&2
		echo "--- See cov_full/log.out for details of what happened ---" 1>&2
	fi
}


# Runs TVC, then copies the .vcf file needed, and deletes the rest of the files generated.
function runTVC {
	mkdir -p ${OUTPUT_DIR}/tvc${TVC_VERSION}_out 2>/dev/null
	checkBamIndex $BAM_FILE ${OUTPUT_DIR}/tvc${TVC_VERSION}_out/log.out

	echo "	$BAM_FILE beginning TVC v${TVC_VERSION} at: `date`"
	#now run TVC v4.2. 
	run_tvc="""${VARIANT_CALLER_DIR}/variant_caller_pipeline.py \
		--input-bam "$BAM_FILE" \
		--reference-fasta "$REF_FASTA" \
		--output-dir "${OUTPUT_DIR}/tvc${TVC_VERSION}_out" \
		--parameters-file "$TVC_JSON" \
		--bin-dir "${VARIANT_CALLER_DIR}" \
		--region-bed  "${TVC_BED}"	"""

	# Add specified options
	if [ "${AMPLISEQ}" == "True" ]; then
		# The --primer-trim-bed will trim the reads to match the amplicons listed in the BED file
		run_tvc="""$run_tvc --primer-trim-bed "${TVC_BED}" """
		# If the ptrim is specified, and it doesn't already exist, then generate it.
		BAM_DIR=`dirname $BAM_FILE`
		if [ "$PTRIM" != "" -a ! "`find ${BAM_DIR}/${PTRIM} -maxdepth 0 -type f 2>/dev/null`" ]; then
			run_tvc="$run_tvc --postprocessed-bam=${BAM_DIR}/${PTRIM} "
		fi
	fi
	if [ "$TVC_HOTSPOT" != "" ]; then
		run_tvc="""$run_tvc --hotspot-vcf "$TVC_HOTSPOT" """
	fi

	# run tvc and write the std_out to the log file
	echo "	running tvc with these options: $run_tvc at `date`" >> ${OUTPUT_DIR}/tvc${TVC_VERSION}_out/log.out
	$run_tvc >> ${OUTPUT_DIR}/tvc${TVC_VERSION}_out/log.out 2>&1
	
	if [ $? -eq 0 ]; then
		#TVC was successful.
		echo "	$BAM_FILE TVC v${TVC_VERSION} finished successfully at: `date`"
	else
		# Something went wrong with TVC. Not copying the data.
		NOERRS="False"
		echo 1>&2
		echo "--- ERROR: $BAM_FILE TVC v${TVC_VERSION} had a problem. Not copying data ---" 1>&2
		echo "--- See tvc${TVC_VERSION}_out/log.out for details of what happened ---" 1>&2
	fi
}


# ----------------------------------------------------------
# ---------------- PROGRAM STARTS HERE ---------------------
# ----------------------------------------------------------


if [ "$RUN_COV" == "True" ]; then
	# Do file checks first
	if [ ! -f $MERGED_BED ]; then
		echo "ERROR: Merged Bed: $MERGED_BED not found. Not running Coverage Analysis."
		NOERRS="False"
	elif [ ! -f $UNMERGED_BED ]; then
		echo "ERROR: Unmerged Bed: $UNMERGED_BED not found. Not running Coverage Anlysis."
		NOERRS="False"
	elif [ "$AMPLISEQ" == "True" -a "`find ${OUTPUT_DIR}/*.amplicon.cov.xls 2>/dev/null`" -a "$FORCED" != "True" ]; then
		echo "	$OUTPUT_DIR already has a .amplicon.cov.xls file. Not running Coverage Analysis. (use --forced to have cov analysis run anyway)"
	elif [ "$TARGETSEQ" == "True" -a "`find ${OUTPUT_DIR}/*.target.cov.xls 2>/dev/null`" -a "$FORCED" != "True" ]; then
		echo "	$OUTPUT_DIR already has a .target.cov.xls file. Not running Coverage Analysis. (use --forced to have cov analysis run anyway)"
	else

		# Only check the flags if TVC or Cov has not yet been run
		if [ "$RUN_TVC" != "True" ]; then
			checkFlags 
		fi
		# If there is no .cov.xls file and the other bed files are found, run coverage analysis
		runCov
	fi
fi
	
if [ "$RUN_TVC" == "True" ]; then
	if [ ! -f $TVC_JSON ]; then 
		echo "ERROR: TVC parameters: $TVC_JSON not found. Not running TVC."
		NOERRS="False"
	elif [ ! -f $TVC_BED ]; then
		echo "ERROR: Project Bed: $TVC_BED not found. Not running TVC."
		NOERRS="False"
	elif [ "$TVC_HOTSPOT" != "" -a ! -f "$TVC_HOTSPOT" ]; then
		echo "ERROR: TVC Hotspot: $TVC_HOTSPOT not found. Not running TVC."
		NOERRS="False"
	elif [ "`find ${OUTPUT_DIR}/${TVC_VERSION}*.vcf 2>/dev/null`" -a "$FORCED" != "True" ]; then
		echo "	$OUTPUT_DIR already has a .vcf file. Not runnning TVC. (use --forced to have tvc run anyway)"
	else
		# Only check the flags if TVC or Cov has not yet been run
		checkFlags 
		# If there is no vcf file present for the current run, and the project BED is found, then run TVC.
		runTVC
	fi
fi

# If there were errors in TVC or Cov analysis, then exit here with an exit status of 1.
if [ "$NOERRS" == "False" ]; then
	echo "Finished with problems."
	exit 1
fi

# If the Cleanup option is specified, cleanup the output
if [ "$CLEANUP" == "True" ]; then
	#echo "Copying the output files and removing the temporary files."
	# Copy the output files we want and remove the temporary files.
	if [ "$TARGETSEQ" == "True" ]; then
		mv ${OUTPUT_DIR}/cov_full/*.target.cov.xls ${OUTPUT_DIR}/ 2>/dev/null
	elif [ "$AMPLISEQ" == "True" ]; then
		mv ${OUTPUT_DIR}/cov_full/*.amplicon.cov.xls ${OUTPUT_DIR}/ 2>/dev/null
	fi
	mv ${OUTPUT_DIR}/tvc${TVC_VERSION}_out/TSVC_variants.vcf ${OUTPUT_DIR}/${TVC_VERSION}_TSVC_variants.vcf 2>/dev/null 

	rm -rf ${OUTPUT_DIR}/cov_full 2>/dev/null
	rm -rf ${OUTPUT_DIR}/tvc${TVC_VERSION}_out 2>/dev/null
#	rm ${BAM_FILE}.bai 2>/dev/null
fi
exit 0
