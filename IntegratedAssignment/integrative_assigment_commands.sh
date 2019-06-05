#set up
export SOFT_DIR=/usr/local/
export WORK_DIR=~/workspace/HTseq/Integrative_Assignment/
export TRIMMOMATIC_JAR=$SOFT_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar
export PICARD_JAR=$SOFT_DIR/picard/picard.jar
export GATK_JAR=$SOFT_DIR/gatk-4.0.1.2/gatk-package-4.0.1.2-local.jar
export GATK_OLD_JAR=~/CourseData/HT_data/software/GenomeAnalysisTK-3.8/GenomeAnalysisTK.jar
export BVATOOLS_JAR=~/CourseData/HT_data/software/bvatools-1.6/bvatools-1.6-full.jar
export REF=$WORK_DIR/reference/


rm -rf $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR
ln -s ~/CourseData/HT_data/Module3/* .

# Parent 1

# Quality
mkdir originalQC_NA12892/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12892/NA12892_CBW_chr1_R2.fastq.gz \
  --threads 2 --regionName ACTL8 --output originalQC_NA12892/


#trim

cat $REF/adapters.fa

mkdir -p reads/NA12892/

java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12892/NA12892_CBW_chr1_R2.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:${REF}/adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12892/NA12892.trim.out

cat reads/NA12892/NA12892.trim.out

mkdir postTrimQC_NA12892/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 2 --regionName ACTL8 --output postTrimQC_NA12892/


# Alignment
mkdir -p alignment/NA12892/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12892\tSM:NA12892\tLB:NA12892\tPU:runNA12892_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/hg19.fa \
  reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx2G -jar ${GATK_JAR} SortSam \
  -I /dev/stdin \
  -O alignment/NA12892/NA12892.sorted.bam \
  -SO coordinate \
  --CREATE_INDEX true \
  --MAX_RECORDS_IN_RAM=500000



# Explore the bam files
samtools view alignment/NA12892/NA12892.sorted.bam | head -n4

# Count the *un-aligned* reads
samtools view -c -f4 alignment/NA12892/NA12892.sorted.bam
# Count the *aligned* reads 
samtools view -c -F4 alignment/NA12892/NA12892.sorted.bam

# Indel realignment
# 1- Find the targets 2- Realign them.
java -Xmx2G  -jar ${GATK_OLD_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/hg19.fa \
  -o alignment/NA12892/realign.intervals \
  -I alignment/NA12892/NA12892.sorted.bam \
  -L chr1

java -Xmx2G -jar ${GATK_OLD_JAR} \
  -T IndelRealigner \
  -R ${REF}/hg19.fa \
  -targetIntervals alignment/NA12892/realign.intervals \
  -o alignment/NA12892/NA12892.realigned.sorted.bam \
  -I alignment/NA12892/NA12892.sorted.bam



# FixMates (optional, see Module 3 main page for explanation)
#java -Xmx2G -jar ${PICARD_JAR} FixMateInformation \
#  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
#  INPUT=alignment/NA12892/NA12892.realigned.sorted.bam \
#  OUTPUT=alignment/NA12892/NA12892.matefixed.sorted.bam

# Mark duplicates
java -Xmx2G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  -I alignment/NA12892/NA12892.realigned.sorted.bam \
  -O alignment/NA12892/NA12892.sorted.dup.bam \
  --METRICS_FILE=alignment/NA12892/NA12892.sorted.dup.metrics

# less alignment/NA12892/NA12892.sorted.dup.metrics


# Recalibration
java -Xmx2G -jar ${GATK_JAR} BaseRecalibrator \
  -R ${REF}/hg19.fa \
  --known-sites ${REF}/dbSNP_135_chr1.vcf.gz \
  -L chr1:17704860-18004860 \
  -O alignment/NA12892/NA12892.sorted.dup.recalibration_report.grp \
  -I alignment/NA12892/NA12892.sorted.dup.bam

java -Xmx2G -jar ${GATK_JAR} ApplyBQSR \
  -R ${REF}/hg19.fa \
  -bqsr alignment/NA12892/NA12892.sorted.dup.recalibration_report.grp \
  -O alignment/NA12892/NA12892.sorted.dup.recal.bam \
  -I alignment/NA12892/NA12892.sorted.dup.bam




##emit 1 warning for the dictionnary
#not working

#java -Xmx2G -jar ${GATK_JAR} \
#  -T PrintReads \
#  -nct 2 \
#  -R ${REF}/hg19.fa \
#  -BQSR alignment/NA12892/NA12892.sorted.dup.recalibration_report.grp \
#  -o alignment/NA12892/NA12892.sorted.dup.recal.bam \
#  -I alignment/NA12892/NA12892.sorted.dup.bam


# Extract Metrics
java  -Xmx2G -jar ${GATK_OLD_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/hg19.fa \
  -o alignment/NA12892/NA12892.sorted.dup.recal.coverage \
  -I alignment/NA12892/NA12892.sorted.dup.recal.bam \
  -L chr1:17700000-18100000
  
## less -S alignment/NA12892/NA12892.sorted.dup.recal.coverage.sample_interval_summary


java -Xmx2G -jar ${GATK_JAR} CollectInsertSizeMetrics \
  -R ${REF}/hg19.fa \
  -I alignment/NA12892/NA12892.sorted.dup.recal.bam \
  -O alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.tsv \
  -H alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.histo.pdf \
  --METRIC_ACCUMULATION_LEVEL LIBRARY



  
## less -S alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.tsv
# head -9  alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.tsv | tail -3 | cut -f5,6

java -Xmx2G -jar ${GATK_JAR} CollectAlignmentSummaryMetrics \
  -R ${REF}/hg19.fa \
  -I alignment/NA12892/NA12892.sorted.dup.recal.bam \
  -O alignment/NA12892/NA12892.sorted.dup.recal.metric.alignment.tsv \
  --METRIC_ACCUMULATION_LEVEL LIBRARY


## less -S alignment/NA12892/NA12892.sorted.dup.recal.metric.alignment.tsv
# head -10  alignment/NA12892/NA12892.sorted.dup.recal.metric.alignment.tsv | tail -4 | cut -f7



  
### parent 2


# Quality
mkdir originalQC_NA12891/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12891/NA12891_CBW_chr1_R2.fastq.gz \
  --threads 2 --regionName ACTL8 --output originalQC_NA12891/


#trim

cat $REF/adapters.fa

mkdir -p reads/NA12891/

java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12891/NA12891_CBW_chr1_R2.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:${REF}/adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12891/NA12891.trim.out

cat reads/NA12891/NA12891.trim.out

mkdir postTrimQC_NA12891/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 2 --regionName ACTL8 --output postTrimQC_NA12891/


# Alignment
mkdir -p alignment/NA12891/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12891\tSM:NA12891\tLB:NA12891\tPU:runNA12891_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/hg19.fa \
  reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx2G -jar ${GATK_JAR} SortSam \
  -I /dev/stdin \
  -O alignment/NA12891/NA12891.sorted.bam \
  -SO coordinate \
  --CREATE_INDEX true \
  --MAX_RECORDS_IN_RAM=500000



# Explore the bam files
samtools view alignment/NA12891/NA12891.sorted.bam | head -n4

# Count the *un-aligned* reads
samtools view -c -f4 alignment/NA12891/NA12891.sorted.bam
# Count the *aligned* reads 
samtools view -c -F4 alignment/NA12891/NA12891.sorted.bam

# Indel realignment
# 1- Find the targets 2- Realign them.
java -Xmx2G  -jar ${GATK_OLD_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/hg19.fa \
  -o alignment/NA12891/realign.intervals \
  -I alignment/NA12891/NA12891.sorted.bam \
  -L chr1

java -Xmx2G -jar ${GATK_OLD_JAR} \
  -T IndelRealigner \
  -R ${REF}/hg19.fa \
  -targetIntervals alignment/NA12891/realign.intervals \
  -o alignment/NA12891/NA12891.realigned.sorted.bam \
  -I alignment/NA12891/NA12891.sorted.bam



# FixMates (optional, see Module 3 main page for explanation)
#java -Xmx2G -jar ${PICARD_JAR} FixMateInformation \
#  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
#  INPUT=alignment/NA12891/NA12891.realigned.sorted.bam \
#  OUTPUT=alignment/NA12891/NA12891.matefixed.sorted.bam

# Mark duplicates
java -Xmx2G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  -I alignment/NA12891/NA12891.realigned.sorted.bam \
  -O alignment/NA12891/NA12891.sorted.dup.bam \
  --METRICS_FILE=alignment/NA12891/NA12891.sorted.dup.metrics

# less alignment/NA12891/NA12891.sorted.dup.metrics


# Recalibration
java -Xmx2G -jar ${GATK_JAR} BaseRecalibrator \
  -R ${REF}/hg19.fa \
  --known-sites ${REF}/dbSNP_135_chr1.vcf.gz \
  -L chr1:17704860-18004860 \
  -O alignment/NA12891/NA12891.sorted.dup.recalibration_report.grp \
  -I alignment/NA12891/NA12891.sorted.dup.bam

java -Xmx2G -jar ${GATK_JAR} ApplyBQSR \
  -R ${REF}/hg19.fa \
  -bqsr alignment/NA12891/NA12891.sorted.dup.recalibration_report.grp \
  -O alignment/NA12891/NA12891.sorted.dup.recal.bam \
  -I alignment/NA12891/NA12891.sorted.dup.bam




##emit 1 warning for the dictionnary
#not working

#java -Xmx2G -jar ${GATK_JAR} \
#  -T PrintReads \
#  -nct 2 \
#  -R ${REF}/hg19.fa \
#  -BQSR alignment/NA12891/NA12891.sorted.dup.recalibration_report.grp \
#  -o alignment/NA12891/NA12891.sorted.dup.recal.bam \
#  -I alignment/NA12891/NA12891.sorted.dup.bam


# Extract Metrics
java  -Xmx2G -jar ${GATK_OLD_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/hg19.fa \
  -o alignment/NA12891/NA12891.sorted.dup.recal.coverage \
  -I alignment/NA12891/NA12891.sorted.dup.recal.bam \
  -L chr1:17700000-18100000
  
## less -S alignment/NA12891/NA12891.sorted.dup.recal.coverage.sample_interval_summary


java -Xmx2G -jar ${GATK_JAR} CollectInsertSizeMetrics \
  -R ${REF}/hg19.fa \
  -I alignment/NA12891/NA12891.sorted.dup.recal.bam \
  -O alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.tsv \
  -H alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.histo.pdf \
  --METRIC_ACCUMULATION_LEVEL LIBRARY



  
## less -S alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.tsv
# head -9  alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.tsv | tail -3 | cut -f5,6

java -Xmx2G -jar ${GATK_JAR} CollectAlignmentSummaryMetrics \
  -R ${REF}/hg19.fa \
  -I alignment/NA12891/NA12891.sorted.dup.recal.bam \
  -O alignment/NA12891/NA12891.sorted.dup.recal.metric.alignment.tsv \
  --METRIC_ACCUMULATION_LEVEL LIBRARY


## less -S alignment/NA12891/NA12891.sorted.dup.recal.metric.alignment.tsv
# head -10  alignment/NA12891/NA12891.sorted.dup.recal.metric.alignment.tsv | tail -4 | cut -f7





