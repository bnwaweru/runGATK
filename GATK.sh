#!/usr/bin/env bash
# Script: GATK.sh
# Usage: bash GATK.sh <Sample_Name> <Reference.fa> <SNPEFF_Db> <KnownDB.vcf>
# Description: bash wrapper for germline variant calling using bwa + GATK
# Does work only with GATK 3.8.1 version
# Contact: s.sivasubramani@cgar.org
# data: 28-02-2018
# set -x

############################
######### TOOLS ############
############################

GATK_PATH="/data/Programs/GATK/GenomeAnalysisTK.jar"
PICARD_PATH="/data/Programs/GATK/picard.jar"
BWA_PATH=$(which bwa)
BCFTOOLS_PATH=$(which bcftools)
SAMTOOLS_PATH=$(which samtools)
BEDTOOLS_PATH=$(which bedtools)
SNPEFF_PATH="/data/Programs/snpEff/snpEff.jar"

#################################
###### Default Parameters #######
#################################

NT=60
NCT=60
# EXP="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0"
EXP="MIN(DP>8) && MIN(MQ>40.00) && MIN(GQ>40.00)"


##############################
###### Default Options #######
##############################

CALCULATEMEM=TRUE
ANNOTATE=FALSE


################################
###### Utility Functions #######
################################


YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'

echo_cmd(){
	echo -e ${YELLOW}CMD: $*${NC} 
	echo "$*"
	eval $*
}

echo_log(){
    text=$*
    size=${#text}
    line=$(printf "%0.s-" $(seq 1 $size))
    echo
    echo -e ${GREEN}$line${NC}
    echo -e ${GREEN}$*${NC}
    echo -e ${GREEN}$line${NC}
    echo
}

isTool(){
	if [ ! -f $* ];
	then
		echo $* does not available.
		exit
	fi
}

isFileExit(){
	if [ ! -f $* ];
	then
		echo $* does not available.
		exit
	fi
}

isFileWarn(){
	if [ ! -f $* ];
	then
		echo $* does not available.
	fi
}

getmem(){
	if [ $CALCULATEMEM == "TRUE" ]
		then
		eval "cat /proc/meminfo | grep MemFree | awk '{print \$2-10000000}'"
	else
		MAXMEM=MAXMEM
	fi
}

echo_usage(){
	echo -e "
	bash GATK.sh -S Sample_Name -s fq -R Reference.fa
	
	Input Files:
	-S Sample name. Script will looks for {Sample_Name}_1.{SUFFIX} and {Sample_Name}_2.{SUFFIX} for mapping. 
	-s Suffix for the read files. fq/fastq/fastq.gz/fq.gz (Tool will look for {Sample_Name}_1.{SUFFIX} and {Sample_Name}_2.{SUFFIX})
	-R Reference FASTA file. .dict file will be created if it does not exist. Path should be having write permission.
	-v Known VCF file for recaliberation. If the given file does not exist, the tool will use called variants for recalliberation.
	
	Annotation:
	-a Enable snpEFF annotation (default:FALSE)
	-d snpEFF database name (No test is enabled as of now. If the mentioned DB does not exist, the tool may try to download. Could be extremly slow based on internet speed)
	
	Variant Calling:
	-e Expression used to filter VCF files (bcftools). (default: \"MIN(DP>8) && MIN(MQ>40.00) && MIN(GQ>40.00)\")
	
	Performance:
	-m Maximum memory for Java VM. Will be over written if -M given.
	-M Let the tool calculate Maximum available memory for JVM (default: TRUE) Not recommended for parallel runs
	-t Number of threads to use for bwa and GATK
	
	FQPN for Tools:
	-G Path for GATK.jar (default: ${GATK_PATH})
	-P Path for Picard.jar (default: ${PICARD_PATH})
	-B Path for bwa (default: $(which bwa || echo "not found, please install bwa"))
	-F Path for bcftools (default: $(which bcftools || echo "not found, please install bcftools"))
	-T Path for samtools (default: $(which samtools || echo "not found, please install samtools"))
	-D Path for snpEFF.jar (default: ${SNPEFF_PATH})
	-E Path for BedTools (default: $(which bedtools || echo "not found, please install bedtools"))

	Warning: 
	1. This script works only with 3.8.1 version of GATK.
	2. Some of the R libraries are needed in making plots (ggplot, gsalib, etc.,).
	3. Make sure snpEff database is present before running the functional annotation.
	4. Make sure the KnownDB.vcf file is sorted accordingly to the reference fasta. (chromosome order)

"
}

##########################
###### Get Options #######
##########################


while getopts ":G:P:B:F:T:E:D:S:s:R:v:d:m:t:e:ahM" option
do
case "${option}"
in
(G)
	GATK_PATH=${OPTARG}
	isTool ${GATK_PATH}
	;;
(P)
	PICARD_PATH=${OPTARG}
	isTool ${PICARD_PATH}
	;;
(B)
	BWA_PATH=${OPTARG}
	isTool ${BWA_PATH}
	;;
(F)
	BCFTOOLS_PATH=${OPTARG}
	isTool ${BCFTOOLS_PATH}
	;;
(T)
	SAMTOOLS_PATH=${OPTARG}
	isTool ${SAMTOOLS_PATH}
	;;
(E)
	BEDTOOLS_PATH=${OPTARG}
	isTool ${BEDTOOLS_PATH}
	;;
(D)
	SNPEFF_PATH=${OPTARG}
	isTool ${SNPEFF_PATH}
	;;
(S)
	SAMPLEID=${OPTARG}
	;;
(s)
	SUFFIX=${OPTARG}
	;;
(R)
	REFERENCE=${OPTARG}
	isFileExit ${REFERENCE} 
	;;
(v)
	KNOWNVCF=${OPTARG}
	isFileWarn ${KNOWNVCF}
	;;
(a)
	ANNOTATE=TRUE
	;;
(d)
	SNPEFFDB=${OPTARG}
	;;
(m)
	MAXMEM=${OPTARG}
	CALCULATEMEM=FALSE
	;;
(M)
	CALCULATEMEM=TRUE
	;;
(t)
	NT=${OPTARG}
	NCT=${OPTARG}
	;;
(e)
	EXP=${OPTARG}
	;;

(h)
	echo_usage
	exit
	;;

(*)
	echo_usage
	exit
	;;
esac
done

if [[ $1 == "" ]];
then
	echo_usage
	exit
fi

if [ ! -f ${REFERENCE%.fa}.bwt ];
then
	echo_log "BWA indexing"
	echo_cmd "${BWA_PATH} index -p ${REFERENCE%.fa} ${REFERENCE}"
fi

echo_log "Performing BWA mem alignment"
echo_cmd "${BWA_PATH} mem -t ${NT} -M -o ${SAMPLEID}.sam ${REFERENCE%.fa} ${SAMPLEID}_1.fq ${SAMPLEID}_2.fq"

if [ ! -f ${REFERENCE%.fa}.dict ];
then
	echo_log "Creating Reference Dictionary.."
	MAXMEM=$(getmem)
	echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${PICARD_PATH} CreateSequenceDictionary R=${REFERENCE} O=${REFERENCE%.fa}.dict"
fi

echo_log "SAM to BAM Convertion + Sort BAM"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${PICARD_PATH} SortSam INPUT=${SAMPLEID}.sam OUTPUT=${SAMPLEID}.sort.bam SORT_ORDER=coordinate"

echo_log "Adding ReadGroups to BAM file"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${PICARD_PATH} AddOrReplaceReadGroups I=${SAMPLEID}.sort.bam O=${SAMPLEID}.rg.sort.bam RGID=${SAMPLEID} RGLB=Pig_ddRAD RGPL=illumina RGPM=NexSeq RGPU=unit1 RGSM=${SAMPLEID}"

echo_log "Collect Alignment Summary"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${PICARD_PATH} CollectAlignmentSummaryMetrics R=${REFERENCE} I=${SAMPLEID}.rg.sort.bam O=${SAMPLEID}.alignment_metrics.txt"

echo_log "Collect Insert Size"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${PICARD_PATH} CollectInsertSizeMetrics INPUT=${SAMPLEID}.rg.sort.bam OUTPUT=${SAMPLEID}.insert_metrics.txt HISTOGRAM_FILE=${SAMPLEID}.insert_size_histogram.pdf"

echo_log "Samtools Depth calculation"
echo_cmd "samtools depth -a ${SAMPLEID}.rg.sort.bam > ${SAMPLEID}.depth_out.txt"

echo_log "Mark Duplicates"
MAXMEM=$(getmem)
echo_cmd "java -jar ${PICARD_PATH} MarkDuplicates INPUT=${SAMPLEID}.rg.sort.bam OUTPUT=${SAMPLEID}.dedup_reads.bam METRICS_FILE=${SAMPLEID}.metrics.txt"

echo_log "Building BAM index"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${PICARD_PATH} BuildBamIndex INPUT=${SAMPLEID}.dedup_reads.bam"

echo_log "Realigner Target Creator"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T RealignerTargetCreator -R ${REFERENCE} -I ${SAMPLEID}.dedup_reads.bam -o ${SAMPLEID}.realignment_targets.list -nt ${NT}"

echo_log "Indel Realigner"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T IndelRealigner -R ${REFERENCE} -I ${SAMPLEID}.dedup_reads.bam -targetIntervals ${SAMPLEID}.realignment_targets.list -o ${SAMPLEID}.realigned_reads.bam"

if [[ ! -f ${KNOWNVCF} ]] || [[ ${KNOWNVCF} == "" ]];
then
	echo "No known variant vcf file is given. Using traditional base recaliberation method by reusing the called SNPs.."
	echo_log "Haplotype caller for caliberation run"
	MAXMEM=$(getmem)
	echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T HaplotypeCaller -R ${REFERENCE} -I ${SAMPLEID}.realigned_reads.bam -o ${SAMPLEID}.raw_variants.vcf -nct ${NCT}"


	echo_log "SNP from Raw VCF"
	MAXMEM=$(getmem)
	echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T SelectVariants -R ${REFERENCE} -V ${SAMPLEID}.raw_variants.vcf -selectType SNP -o ${SAMPLEID}.raw_snps.vcf"

	echo_log "InDels from Raw VCF"
	MAXMEM=$(getmem)
 	echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T SelectVariants -R ${REFERENCE} -V ${SAMPLEID}.raw_variants.vcf -selectType INDEL -o ${SAMPLEID}.raw_indels.vcf"

	echo_log "Filtering SNP VCF 1"
	MAXMEM=$(getmem)
#	echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T VariantFiltration -R ${REFERENCE} -V ${SAMPLEID}.raw_snps.vcf --filterExpression \"${EXP}\" --filterName \"basic_snp_filter\" -o ${SAMPLEID}.filtered_snps.vcf"
	echo_cmd "${BCFTOOLS_PATH} view -i 'MIN(DP>8) && MIN(MQ>40.00) && MIN(GQ>40.00)' ${SAMPLEID}.raw_snps.vcf -o ${SAMPLEID}.filtered_snps.vcf"

	echo_log "Filtering InDel VCF"
	MAXMEM=$(getmem)
#	echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T VariantFiltration -R ${REFERENCE} -V ${SAMPLEID}.raw_indels.vcf --filterExpression \"${EXP}\" --filterName \"basic_indel_filter\" -o ${SAMPLEID}.filtered_indels.vcf"
	echo_cmd "${BCFTOOLS_PATH} view -i 'MIN(DP>8) && MIN(MQ>40.00) && MIN(GQ>40.00)' ${SAMPLEID}.raw_indels.vcf -o ${SAMPLEID}.filtered_indels.vcf"

fi

if [[ -f ${KNOWNVCF} ]] && [[ ${KNOWNVCF} != "" ]];
then
	echo_log "SNPs from reference VCF"
	MAXMEM=$(getmem)
	echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T SelectVariants -R ${REFERENCE} -V ${KNOWNVCF} -selectType SNP -o ${SAMPLEID}.filtered_snps.vcf"

	echo_log "InDel from reference VCF"
	MAXMEM=$(getmem)
	echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T SelectVariants -R ${REFERENCE} -V ${KNOWNVCF} -selectType INDEL -o ${SAMPLEID}.filtered_indels.vcf"
fi

echo_log "Base recaliberation 1"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T BaseRecalibrator -R ${REFERENCE} -I ${SAMPLEID}.realigned_reads.bam -knownSites ${SAMPLEID}.filtered_snps.vcf -knownSites ${SAMPLEID}.filtered_indels.vcf -o ${SAMPLEID}.recal_data.table -nct ${NCT}"

echo_log "Base recaliberation 2"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T BaseRecalibrator -R ${REFERENCE} -I ${SAMPLEID}.realigned_reads.bam -knownSites ${SAMPLEID}.filtered_snps.vcf -knownSites ${SAMPLEID}.filtered_indels.vcf -BQSR ${SAMPLEID}.recal_data.table -o post_${SAMPLEID}.recal_data.table -nct ${NCT}"

echo_log "Analyze covariates"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T AnalyzeCovariates -R ${REFERENCE} -before ${SAMPLEID}.recal_data.table -after post_${SAMPLEID}.recal_data.table -plots ${SAMPLEID}.recalibration_plots.pdf"

echo_log "Print reads"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T PrintReads -R ${REFERENCE} -I ${SAMPLEID}.realigned_reads.bam -BQSR ${SAMPLEID}.recal_data.table -o ${SAMPLEID}.recal_reads.bam -nct ${NCT}"

echo_log "Variant Calling (HaplotypeCaller)"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T HaplotypeCaller -R ${REFERENCE} -I ${SAMPLEID}.recal_reads.bam -o ${SAMPLEID}.raw_variants_recal.vcf  -nct ${NCT}"

echo_log "SNP from Raw VCF (recalibrated)"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T SelectVariants -R ${REFERENCE} -V ${SAMPLEID}.raw_variants_recal.vcf -selectType SNP -o ${SAMPLEID}.raw_snps_recal.vcf"

echo_log "InDel from Raw VCF (recalibrated)"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T SelectVariants -R ${REFERENCE} -V ${SAMPLEID}.raw_variants_recal.vcf -selectType INDEL -o ${SAMPLEID}.raw_indels_recal.vcf"

echo_log "Filtering SNP VCF (recalibrated)"
MAXMEM=$(getmem)
# echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T VariantFiltration -R ${REFERENCE} -V ${SAMPLEID}.raw_snps_recal.vcf --filterExpression \"${EXP}\" --filterName \"basic_snp_filter\" -o ${SAMPLEID}.filtered_snps_final.vcf"
echo_cmd "${BCFTOOLS_PATH} view -i 'MIN(DP>8) && MIN(MQ>40.00) && MIN(GQ>40.00)' ${SAMPLEID}.raw_snps_recal.vcf -o ${SAMPLEID}.filtered_snps_final.vcf"


echo_log "Filtering InDel VCF (recalibrated)"
MAXMEM=$(getmem)
# echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${GATK_PATH} -T VariantFiltration -R ${REFERENCE} -V ${SAMPLEID}.raw_indels_recal.vcf --filterExpression \"${EXP}\" --filterName \"basic_indel_filter\" -o ${SAMPLEID}.filtered_indels_recal.vcf"
echo_cmd "${BCFTOOLS_PATH} view -i 'MIN(DP>8) && MIN(MQ>40.00) && MIN(GQ>40.00)' ${SAMPLEID}.raw_indel_recal.vcf -o ${SAMPLEID}.filtered_indels_recal.vcf"


echo_log "Functional Annotation using SnpEff"
MAXMEM=$(getmem)
echo_cmd "java -Xmx${MAXMEM%??????}g -jar ${SNPEFF_PATH} eff -htmlStats -csvStats ${SAMPLEID}.snpEff.stats.csv -v ${SNPEFFDB} ${SAMPLEID}.filtered_snps_final.vcf > ${SAMPLEID}.filtered_snps_final.ann.vcf"
echo_cmd "mv snpEff_summary.html ${SAMPLEID}.snpEff.stats.html"

echo_log "BedGraph construction"
echo_cmd "${BEDTOOLS_PATH} genomecov -bga -ibam ${SAMPLEID}.recal_reads.bam > ${SAMPLEID}.genomecov.bedgraph"

