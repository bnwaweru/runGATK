# runGATK
bash wrapper for running GATK (3.8.1) pipeline

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
	-e Expression used to filter VCF files (bcftools). (default: "MIN(DP>8) && MIN(MQ>40.00) && MIN(GQ>40.00)")
	
	Performance:
	-m Maximum memory for Java VM. Will be over written if -M given.
	-M Let the tool calculate Maximum available memory for JVM (default: TRUE) Not recommended for parallel runs
	-t Number of threads to use for bwa and GATK
	
	FQPN for Tools:
	-G Path for GATK.jar (default: /data/Programs/GATK/GenomeAnalysisTK.jar)
	-P Path for Picard.jar (default: /data/Programs/GATK/picard.jar)
	-B Path for bwa (default: /opt/anaconda2/bin/bwa)
	-F Path for bcftools (default: /opt/anaconda2/bin/bcftools)
	-T Path for samtools (default: /opt/anaconda2/bin/samtools)
	-D Path for snpEFF.jar (default: /data/Programs/snpEff/snpEff.jar)
	-E Path for BedTools (default: not found, please install bedtools)

	Warning: 
	1. This script works only with 3.8.1 version of GATK.
	2. Some of the R libraries are needed in making plots (ggplot, gsalib, etc.,).
	3. Make sure snpEff database is present before running the functional annotation.
	4. Make sure the KnownDB.vcf file is sorted accordingly to the reference fasta. (chromosome order)


