#! /bin/bash

# GOAL: change the flags of runs in a bam file.

# Parameters for Qsub:
#$ -S /bin/bash
#$ -N flags 
#$ -o flags.log
#$ -e flags.log
#$ -cwd
#$ -V

PICARD_DIR="/opt/picard/picard-tools-current"

function usage {
cat << EOF
USAGE: bash flags.sh <path/to/BAM_File.bam> [OPTIONS]
Options:
	-r | --remove_dup	(Default action. Will use samtools to check if there are duplicate flags. If there are, remove the duplicate flags from the bam file using Picard tools)
	-f | --flag_dups	(Will place the duplicate flag in the bam file)
	-cl | --cleanup		(Will replace the original BAM file with the modified bam. Will also delete the bam's index file)	
EOF
exit 8
}

# Default parameters
OPTION="REMOVE_DUP_FLAGS"
RUNNING="Running flags.sh at `date` with these options: "

if [ $# -eq 0 ]; then
	usage
fi

while :
do
	case $1 in
		*.bam)
			BAM_FILE=$1
			RUNNING="$RUNNING BAM: $1, "
			shift 
			;;
		-r | --remove_dup)
			OPTION="REMOVE_DUP_FLAGS"
			shift
			;;
		-f | --flag_dups)
			OPTION="FLAG_DUPS"
			RUNNING="$RUNNING --flag_dups, "
			shift
			;;
		-cl | --cleanup)
			CLEANUP="True"
			RUNNING="$RUNNING --cleanup, "
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
if [ "$OPTION" == "REMOVE_DUP_FLAGS" ]; then
	RUNNING="$RUNNING --remove_dup"
fi
echo "	$RUNNING"

# Now check the files to see if the bam file exists.
if [ "$BAM_FILE" == '' ]; then
	echo "USAGE: MUST SPECIFY A BAM FILE" 1>&2
#	usage
	exit 4
elif [ ! "`find $BAM_FILE -maxdepth 0 2>/dev/null`" ]; then
	echo "ERROR: BAM file: $BAM_FILE not found" 1>&2
	exit 4
fi

# --------------- PROGRAM STARTS HERE ---------------------


# first use flagstat to tell us if the bam file has flags or not.
#flagstat=`samtools flagstat $BAM_FILE | head -n 2 | tail -n 1 | grep -oE "^[0-9]+"`
flagstat=`samtools view -f 0x0400 $BAM_FILE | wc -l`

if [ "$OPTION" == "REMOVE_DUP_FLAGS" ]; then
	# Check to see if there are duplicate flags in this bam file. If there are, then remove them.
	if [ $flagstat -ne 0 ]; then
		echo "	$BAM_FILE Removing duplicate flags at: `date`."
		bam_file_name=`basename $BAM_FILE`
		bam_dir=`dirname $BAM_FILE`
		# Run picard tools to remove the flags
		java -jar ${PICARD_DIR}/RevertSam.jar \
			INPUT=$BAM_FILE \
			OUTPUT=${bam_dir}/NoFlags_${bam_file_name} \
			REMOVE_ALIGNMENT_INFORMATION=FALSE \
			REMOVE_DUPLICATE_INFORMATION=TRUE \
			RESTORE_ORIGINAL_QUALITIES=FALSE

		MODIFIED_BAM="${bam_dir}/NoFlags_${bam_file_name}"
	else
		echo "	$BAM_FILE has no duplicate flags."
	fi
elif [ "$OPTION" == "FLAG_DUPS" ]; then
	if [ $flagstat -eq 0 ]; then
		bam_file_name=`basename $BAM_FILE`
		bam_dir=`dirname $BAM_FILE`
		echo "	$BAM_FILE Marking Duplicates at: `date`"
		# MarkDuplicates marks the reads that are duplicates as duplicates in the bam file.
		# METRICS_FILE is the file to write duplication metrics to. It is required, but we can delete it after
		java -jar ${PICARD_DIR}/MarkDuplicates.jar \
			INPUT=${BAM_FILE} \
			OUTPUT=${bam_dir}/Flag_${bam_file_name} \
			METRICS_FILE=metrics.txt

		MODIFIED_BAM="${bam_dir}/Flag_${bam_file_name}"
	else
		echo "	$BAM_FILE already has $flagstat duplicate flags."
	fi
fi


# If cleanup is specified, the new bam file will replace the original bam file.
if [ "$CLEANUP" == "True" -a "$MODIFIED_BAM" != "" ]; then
	# The bam file probably needs to be reindexed after making these changes.
	rm ${BAM_FILE}.bai 2>/dev/null
	mv $MODIFIED_BAM $BAM_FILE
	rm metrics.txt 2>/dev/null
fi
echo "	Finished at: `date`"

