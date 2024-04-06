#!/usr/bin/bash

usage() {

	echo "Usage: $0 -a <ascession number> -d <Directory_To_Place_FastQ_Files> -r <reference_fasta> -s <sample_name> [<GATK Wrapper Script Directory>]"
	exit 1
}

#parsing the options

while getopts "a:d:r:s:g:" opt; do 
	case $opt in 
		a) 
			SRA_ACC="$OPTARG"
			;;
		d) 
			FASTQ_DIR="$OPTARG"
			;;
		r)
			REF="$OPTARG"
			;;
		s)
			SAMPLENAME="$OPTARG"
			;;
		g)
			GATK_DIR="$OPTARG"
			;;
		\?)
			echo "Invalid Option: -$OPTARG" >&2
			usage
			;;
		:)
			echo "Option -$OPTARG requires argument." >&2
			usage
			;;
	esac
done

if [ -z "$SRA_ACC" ] || [ -z "$FASTQ_DIR" ] ||  [ -z "$REF" ] || [ -z "$SAMPLENAME" ]; then
	echo "All options (-a -d -r -s) are required."
	usage
fi

if [ -z "$GATK_DIR"]; then
	GATK_DIR="."
fi

if ! command -v fastq-dump &> /dev/null; then 
	echo "Error: fastq-dump command not found. Make sure SRATools is installed." >&2
	exit 1
fi

if [ ! -d "$FASTQ_DIR" ]; then 
	echo "Directory you have defined not there. Creating directory now..."
	mkdir -p "$FASTQ_DIR"
fi

echo "Running fastq-dump to download and convert SRA to fastq..."

fastq-dump --split-3 $SRA_ACC

if [ $? -ne 0]; then 
	echo "Error fastqq-dump has not worked properly" >&2
	exit 1
fi

mv ${SRA_ACC}* $FASTQ_DIR

cd $FASTQ_DIR || exit

echo "Running GATK Wrapper Script..."
"$GATK_DIR/GATK_no_scattergather.sh" -r $REF -b $FASTQ_DIR -s $SAMPLE_NAME 

if [ $? -ne 0]; then
	echo "Error: GATK wrapper script failed." >&2
	exit 1
fi

echo "$0 has finished.."

			
