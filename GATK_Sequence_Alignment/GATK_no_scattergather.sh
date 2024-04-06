#!/usr/bin/bash 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source $DIR/_Bash_functions/File_input_output.sh

prepareuBams() {
	localDIR=$1
	SAMPLE=$2
	GATK=$3
	PAIRED=$4

	# FASTQ files (Single-end or Paired-end)
	R1=
	R2=
	if [ $PAIRED == 'y' ] ; then
		findFastqPairsInFolder $localDIR
		[[ -f $R1 ]] || die "Haven't found R1 in $localDIR"
		[[ -f $R2 ]] || die "Haven't found R2 in $localDIR"
	else
		findFastqInFolder $localDIR
		[[ -f $R1 ]] || die "Haven't found R1 in $localDIR"
	fi

	# Output
	OUTPUT="${localDIR}/${SAMPLE}.ubam.bam"
	#echo $localDIR
	#echo $SAMPLE

	# Paired
	if [ $PAIRED == 'y' ] ; then
		$GATK FastqToSam \
		--FASTQ $R1 \
		--FASTQ2 $R2 \
		--OUTPUT $OUTPUT \
		--READ_GROUP_NAME 1 \
		--SAMPLE_NAME $SAMPLE \
		--LIBRARY_NAME NextSeq_Illumina \
		--PLATFORM illumina
	else
		$GATK FastqToSam \
		--FASTQ $R1 \
		--OUTPUT $OUTPUT \
		--READ_GROUP_NAME 1 \
		--SAMPLE_NAME $SAMPLE \
		--LIBRARY_NAME NextSeq_Illumina \
		--PLATFORM illumina
	fi
}

preprocessing() {
	echo "Initiating data preprocessing..."
	ref_fasta=$1
	unmapped_bam=$2
	gatk_path=$3
	picard=$4
	PAIRED=$5

	# Input/Output
	#unmapped_bam="${localDIR}/${SAMPLE}.ubam.bam"
	SamToFastqAndBwaMem_output_bam_basename="${unmapped_bam}.unmerged"
	MergeBamAlignment_output_bam_basename="${unmapped_bam}.aligned.unsorted"
	MarkDuplicates_output_bam_basename="${unmapped_bam}.aligned.unsorted.duplicates_marked"
	SortAndFixTags_output_bam_basename="${unmapped_bam}.aligned.duplicate_marked.sorted"
	metrics_filename="${unmapped_bam}.duplicate_metrics"

	# bwa command
	bwa_commandline=
	if [ $PAIRED == 'y' ] ; then
		bwa_commandline="bwa mem -M -p -v 3 -t 16 $ref_fasta"
	else
		bwa_commandline="bwa mem -M -v 3 -t 16 $ref_fasta"
	fi

	# BWA version
	bwa_version=`bwa 2>&1 | \
	grep -e '^Version' | \
	sed 's/Version: //'`

	# SamToFastqAndBwaMem
	${picard} SamToFastq \
	-INPUT ${unmapped_bam} \
	-FASTQ /dev/stdout \
	-INTERLEAVE true \
	-NON_PF true \
	| \
	${bwa_commandline} /dev/stdin 2> >(tee ${unmapped_bam}.bwa.stderr.log >&2) > ${SamToFastqAndBwaMem_output_bam_basename}.sam
	samtools view -1S ${SamToFastqAndBwaMem_output_bam_basename}.sam > ${SamToFastqAndBwaMem_output_bam_basename}.bam

	# MergeBamAlignment
	${gatk_path} MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ALIGNED_BAM ${SamToFastqAndBwaMem_output_bam_basename}.bam \
    --UNMAPPED_BAM ${unmapped_bam} \
    --OUTPUT ${MergeBamAlignment_output_bam_basename}.bam \
    --REFERENCE_SEQUENCE ${ref_fasta} \
    --PAIRED_RUN true \
    --SORT_ORDER "unsorted" \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --CLIP_ADAPTERS false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --PROGRAM_RECORD_ID "bwamem" \
    --PROGRAM_GROUP_VERSION "${bwa_version}" \
    --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
    --PROGRAM_GROUP_NAME "bwamem" \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true

    # MarkDuplicates
    ${gatk_path} MarkDuplicates \
	--INPUT ${MergeBamAlignment_output_bam_basename}.bam \
	--OUTPUT ${MarkDuplicates_output_bam_basename}.bam \
	--METRICS_FILE ${metrics_filename} \
	--VALIDATION_STRINGENCY SILENT \
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
	--ASSUME_SORT_ORDER "queryname" \
	--CREATE_MD5_FILE true

	# SortAndFixTags (Sort aggregated+deduped BAM file and fix tags)
	${gatk_path} SortSam \
	--INPUT ${MarkDuplicates_output_bam_basename}.bam \
	--OUTPUT /dev/stdout \
	--SORT_ORDER "coordinate" \
	--CREATE_INDEX false \
	--CREATE_MD5_FILE false \
	| \
	#${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_fix}" \
	${gatk_path} \
	SetNmMdAndUqTags \
	--INPUT /dev/stdin \
	--OUTPUT ${SortAndFixTags_output_bam_basename}.bam \
	--CREATE_INDEX true \
	--CREATE_MD5_FILE true \
	--REFERENCE_SEQUENCE ${ref_fasta}

	#File analysis_ready_bam = SortAndFixTags.output_bam
	#File analysis_ready_bam_index = SortAndFixTags.output_bam_index

	echo "Preprocessing completed."
}

generateIntervals() {
	REF=$1
	localDIR=$2
	GATK=$3

	# Input files
	DICT=${REF%.*}.dict

	echo "Generating intervals for parallel analysis..."

	awk '{print $1"\t1\t"$2"\t+\t"$1}' $REF.fai | cat $DICT - > $localDIR/hc.interval_list
	awk '{print $1":1-"$2}' $REF.fai > $localDIR/postprocessing.intervals

	mkdir $localDIR/scattered_intervals

	#$GATK --java-options "-Xmx10G" IntervalListTools \
	#$GATK -Xmx10G IntervalListTools \
    $GATK IntervalListTools \
	--INPUT $localDIR/hc.interval_list \
    --OUTPUT $localDIR/scattered_intervals \
    --UNIQUE true \
    --SORT true \
    --SCATTER_COUNT 3 \
    --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --PADDING 0 \
    --ACTION CONCAT \
    --INCLUDE_FILTERED false \
    --INVERT false \
    --VERBOSITY INFO \
    --QUIET false \
    --VALIDATION_STRINGENCY STRICT \
    --COMPRESSION_LEVEL 5 \
    --MAX_RECORDS_IN_RAM 500000 \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false

	ls $localDIR/scattered_intervals/*/* > $localDIR/scatter_list.txt
	echo "Intervals successfully created."
}

haplotypecaller() {
	echo "Initiating haplotype calling..."
	ref_fasta=$1
	localDIR=$2
	SAMPLE=$3
	ploidy_to_use=$4
	gatk_path=$5

	# Input/Output files
	input_bam="${localDIR}/${SAMPLE}.ubam.bam.aligned.duplicate_marked.sorted.bam" 
	#DICT=${ref_fasta%.*}.dict
	#BAI=${input_bam%.*}.bai
	scatter_list="${localDIR}/scatter_list.txt"
	output_filename="${localDIR}/${SAMPLE}."
	output_suffix=".g.vcf.gz"

	# HaplotypeCaller
	i=1
	output_VCF_list=""
	while IFS="" read -r interval_list || [ -n "$p" ]
	do
		#echo $interval_list
		output_vcf="${output_filename}${i}${output_suffix}"
		output_VCF_list="${output_VCF_list} --INPUT ${output_vcf}"
		#echo "output_vcf_list = ${output_VCF_list}"

		${gatk_path} \
		HaplotypeCaller \
		-R ${ref_fasta} \
		-I ${input_bam} \
		-L ${interval_list} \
		-O "${output_vcf}" \
		-ploidy ${ploidy_to_use} \
		-contamination 0 \
		-ERC GVCF

		i=$((i+1))
	done < ${scatter_list}

	# Merge VCFs
	output_vcf="${output_filename}merged${output_suffix}"

	${gatk_path} \
	MergeVcfs \
    ${output_VCF_list} \
    --OUTPUT ${output_vcf}

	echo "Haplotype calling completed."
}

postprocessing() {
	echo "Initiating data postprocessing..."
	ref_fasta=$1
	localDIR=$2
	SAMPLE=$3
	gatk_path=$4

	# Input/Output
	unpadded_intervals_file="${localDIR}/postprocessing.intervals"
	input_gvcfs="${localDIR}/${SAMPLE}.merged.g.vcf.gz"

	# Process intervals
	num_of_original_intervals=0
	num_gvcfs=1
	while IFS="" read -r interval_list || [ -n "$p" ]
	do
		#echo $interval_list
		num_of_original_intervals=$((num_of_original_intervals+1))
	done < ${unpadded_intervals_file}

	# Make a 2.5:1 interval number to samples in callset ratio interval list
	# bc (basic calculator) needed for float point calculations. The option -l is equivalent to --mathlib; it loads the standard math library.
	# Enclosing the whole expression between double parenthesis (( )) will translate these values to respectively true or false.
	possible_merge_count=$(echo "scale=2;$num_of_original_intervals/$num_gvcfs/2.5" | bc)
	merge_count=$possible_merge_count
	if ((  $(echo "$possible_merge_count < 1" | bc -l ) )) ; then
		merge_count=1
	fi
	echo "merge count = ${merge_count}. Making new intervals file using python..."

	# Make new intervals file using python (not sure this does anything!)
	intervals=${unpadded_intervals_file}
	output_intervals="${unpadded_intervals_file}.parsed"
	python << END
def parse_interval(interval):
    colon_split = interval.split(":")
    chromosome = colon_split[0]
    dash_split = colon_split[1].split("-")
    start = int(dash_split[0])
    end = int(dash_split[1])
    return chromosome, start, end

def add_interval(chr, start, end):
        lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        return chr, start, end

count = 0
chain_count = ${merge_count}
l_chr, l_start, l_end = "", 0, 0
lines_to_write = []
with open("${intervals}") as f:
    with open("${output_intervals}", "w") as f1:
        for line in f.readlines():
            # initialization
            if count == 0:
                w_chr, w_start, w_end = parse_interval(line)
                count = 1
                continue
            # reached number to combine, so spit out and start over
            if count == chain_count:
                l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                w_chr, w_start, w_end = parse_interval(line)
                count = 1
                continue

            c_chr, c_start, c_end = parse_interval(line)
            # if adjacent keep the chain going
            if c_chr == w_chr and c_start == w_end + 1:
                w_end = c_end
                count += 1
                continue
            # not adjacent, end here and start a new chain
            else:
                l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                w_chr, w_start, w_end = parse_interval(line)
                count = 1
        if l_char != w_chr or l_start != w_start or l_end != w_end:
            add_interval(w_chr, w_start, w_end)
        f1.writelines("\n".join(lines_to_write))
END

	# for each interval
	batch_size=5
	workspace_dir_name="${localDIR}/genomicsdb"
	input_list="${localDIR}/inputs.list"
	GenotypeGVCFs_output_vcf_filename="${localDIR}/output.vcf.gz"
	echo -e "${SAMPLE}\t${input_gvcfs}" > ${input_list}
	mkdir "${localDIR}/scattered_vcfs"
	output_VCF_list=""
	while IFS="" read -r interval || [ -n "$p" ]
	do
		echo "postprocessing ${interval}"

		# ImportGVCFs
		CMD1="${gatk_path} \
    	GenomicsDBImport \
    	--genomicsdb-workspace-path ${workspace_dir_name} \
    	--batch-size ${batch_size} \
    	-L ${interval} \
    	--sample-name-map ${input_list} \
    	--reader-threads 5 \
    	-ip 500"
    	echo "CMD1: $CMD1"
    	eval $CMD1

    	# tar and untar?
    	CMD2="tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}"
    	echo "CMD2: $CMD2"
    	eval $CMD2

    	workspace_tar=${workspace_dir_name}.tar
    	CMD3="tar -xf ${workspace_tar}"
    	echo "CMD3: $CMD3"
    	eval $CMD3

    	CMD4="${gatk_path} \
    	GenotypeGVCFs \
    	-R ${ref_fasta} \
    	-O ${GenotypeGVCFs_output_vcf_filename} \
    	-G StandardAnnotation \
    	--only-output-calls-starting-in-intervals \
    	--use-new-qual-calculator \
    	--include-non-variant-sites \
    	-V gendb://${workspace_dir_name} \
    	-L ${interval}"
    	echo "CMD4: $CMD4"
    	eval $CMD4

    	# GenotypeGVCFs
    	variant_filtered_vcf_filename="${localDIR}/scattered_vcfs/$SAMPLE.${interval}.variant_filtered.vcf.gz"
        sites_only_vcf_filename="${localDIR}/scattered_vcfs/$SAMPLE.${interval}.sites_only.variant_filtered.vcf.gz"
    	CMD5="${gatk_path} \
    	VariantFiltration \
    	-filter \"QD < 2.0\" --filter-name \"QD\" \
    	-filter \"FS > 60.0\" --filter-name \"FS\" \
    	-filter \"MQ < 40.0\" --filter-name \"MQ\" \
    	-O ${variant_filtered_vcf_filename} \
    	-V ${GenotypeGVCFs_output_vcf_filename}"
    	echo "CMD5: $CMD5"
    	eval $CMD5

    	output_VCF_list="${output_VCF_list} --input ${variant_filtered_vcf_filename}"

    	CMD6="${gatk_path} \
    	MakeSitesOnlyVcf \
    	--INPUT ${variant_filtered_vcf_filename} \
    	--OUTPUT ${sites_only_vcf_filename}"
    	echo "CMD6: $CMD6"
    	eval $CMD6

    	# Clean up
    	CMD7="rm -r ${workspace_dir_name}"
    	echo "CMD7: $CMD7"
    	eval $CMD7

	done < ${output_intervals}

	# GatherVcfs as SitesOnlyGatherVcf
	GatherVcfsCloud_output_vcf_filename="${localDIR}/${SAMPLE}.vcf.gz"
	CMD8="${gatk_path} \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    ${output_VCF_list} \
    --output ${GatherVcfsCloud_output_vcf_filename}"
    echo "CMD8: $CMD8"
    eval $CMD8

    CMD9="${gatk_path} \
    IndexFeatureFile \
    --feature-file ${GatherVcfsCloud_output_vcf_filename}"
    echo "CMD9: $CMD9"
    eval $CMD9
	
}

qualityCheck() {
	echo "Quality check begins..."
	localDIR=$1
	SAMPLE=$2

	# Input files
	VCF="${localDIR}/${SAMPLE}.vcf.gz"
	echo $VCF

	if [[ ! -r $VCF ]]; then
		warn "QualityCheck: The following pipeline finished unsuccessfully: $VCF"
	fi

	gunzip -c $VCF > ${localDIR}/${SAMPLE}.vcf
	grep -E '(PASS|^#)' ${localDIR}/${SAMPLE}.vcf > ${localDIR}/${SAMPLE}_filtered.vcf
	#rm ${localDIR}/${SAMPLE}.vcf
	gzip ${localDIR}/${SAMPLE}_filtered.vcf

	echo "Quality check finished."
}

diemsg() { 
	echo "Usage: $0 -r <reference.fasta> -b <folder with paired FASTQ (full directory required)> -s <sample name>"
	echo "Optional: -p Ploidy (1/2/3) [1]"
	echo "          -i ISCA project ID. To be implemented"
	echo "          -d Paired reads (y/n) [y]"
	echo ""
	echo "Custom (GATK v4.4.4Java-1.8+, BWA, and PicardTools paths can be defined at the beginning of main.sh."
	echo ""
	echo "Notes: Requires >8G memory. For local processing docker is required."
	die ""
}
[ $# -gt 1 ] || diemsg

# Gather command line settings in opts (a, b, c, etc.)
declare -A opts
declare -A BASH_EXE
declare -A WDL_PATHS
opts[p]="1"
opts[d]="y"
getOpts $@

# Check input file(s) are readable
[[ -r ${opts[r]} ]] || die "File: ${opts[r]} does not appear to be valid"
[[ -d ${opts[b]} ]] || die "Folder: ${opts[b]} does not appear to be valid"
[[ -n ${opts[s]} ]] || die "Undefined sample name" 
[[ ${opts[p]} =~ ^(1|2|3)$ ]] || die "Ploidy must equal 1 or 2 or 3 (no text)"
[[ ${opts[d]} =~ ^(y|n)$ ]] || die "Invalid opt_d settings specified (case sensitive)"

echo "Running $0 -r ${opts[r]} -b ${opts[b]}  -s ${opts[s]} -p ${opts[p]} -d ${opts[d]}"

# Paths to programs (used to be server specific including maxwell)
BASH_EXE=(
	["gatk"]="java -Dsamjdk.compression_level=2 -jar $DIR/../java/gatk-package-4.4.0.0-local.jar"
	["bwa"]="bwa"
	["samtools"]="samtools"
)

WDL_PATHS=(
	["gitc"]="null"
	["gatk"]="java -Dsamjdk.compression_level=2 -jar $DIR/../java/gatk-package-4.4.0.0-local.jar"
	["picard"]="java -Dsamjdk.compression_level=5 -Dpicard.useLegacyParser=false -Xms3000m -jar $DIR/../java/picard.jar"
)

# Index reference 1
if [ ! -r ${opts[r]}.amb ]; then
	echo "BWA: Indexing reference and remove any previously associated file..."
	${BASH_EXE[bwa]} index -a is ${opts[r]}
fi

# Index reference 2
if [ ! -r ${opts[r]}.fai ]; then
	echo "SAMTOOLS: Creating reference index file..."
	${BASH_EXE[samtools]} faidx ${opts[r]}
fi

# Create FASTA dictionary
DICT=${opts[r]%.*}.dict
if [ ! -r $DICT ]; then
	echo "GATK: Creating reference dictionary..."
	${BASH_EXE[gatk]} CreateSequenceDictionary -R ${opts[r]}
fi

echo "Preparing ubams..."
prepareuBams ${opts[b]} ${opts[s]} "${BASH_EXE[gatk]}" ${opts[d]}

echo "Preprocessing ubams..."
preprocessing ${opts[r]} "${opts[b]}/${opts[s]}.ubam.bam" "${WDL_PATHS[gatk]}" "${WDL_PATHS[picard]}" ${opts[d]}

echo "Generate Intervals..."
generateIntervals ${opts[r]} ${opts[b]} "${BASH_EXE[gatk]}"

echo "HaplotypeCaller..."
haplotypecaller ${opts[r]} ${opts[b]} ${opts[s]} ${opts[p]} "${WDL_PATHS[gatk]}"

echo "Postprocessing..."
postprocessing ${opts[r]} ${opts[b]} ${opts[s]} "${WDL_PATHS[gatk]}"

echo "Final Actions..."
qualityCheck ${opts[b]} ${opts[s]}

