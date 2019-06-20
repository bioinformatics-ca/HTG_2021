---
layout: tutorial_page
permalink: /htseq_2019_IA
title: HTseq Integrated Assignment
header1: Workshop Pages for Students
header2: Informatics for High-throughput Sequencing Data Analysis 2019 Integrated Assignment
image: /site_images/CBW_High-throughput_icon.jpg
home: https://bioinformaticsdotca.github.io/htseq_2019
---


-----------------------

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

-----------------------

# CBW HT-seq Integrative Assignment

 
Written originally by Mathieu Bourgey, edited by Florence Cavalli


### Task
We will perform the same analysis as in Module 3 but using the mother and father samples i.e sample NA12891 and NA12891.

```bash
The fastq files are in the following directory of the cloud instance: ~/CourseData/HT_data/Module3/

 * raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz
 * raw_reads/NA12891/NA12891_CBW_chr1_R2.fastq.gz
 * raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz
 * raw_reads/NA12892/NA12892_CBW_chr1_R2.fastq.gz
```


### Environment setup

```bash
# Start a fresh docker container:

docker run --privileged -v /tmp:/tmp --network host -it -w $PWD -v $HOME:$HOME -v /media:/media --user $UID:$GROUPS -v /etc/group:/etc/group -v /etc/passwd:/etc/passwd c3genomics/genpipes:0.8

```

Then inside the docker container:

```bash
#set up environment variables
export REF=$WORK_DIR/reference/
export WORK_DIR=~/workspace/HTseq/Integrative_Assignment/

rm -rf $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR
ln -fs ~/CourseData/HT_data/Module3/* .

# Load the software modules

module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/trimmomatic/0.36 mugqic/samtools/1.9 mugqic/bwa/0.7.17 mugqic/GenomeAnalysisTK/4.1.0.0 mugqic/R_Bioconductor/3.5.0_3.7

```


Task list:

1. Check read QC

2. Trim unreliable bases from the read ends

3. Align the reads to the reference

4. Sort the alignments by chromosome position

5. Realign short indels

6. Fixe mate issues (optional)

7. Mark duplicates

8. Recalibrate the Base Quality

9. Generate alignment metrics


Discussion/Questions:

1. Explain the purpose of each step

2. Which software tool can be used for each step 




The full commands can be downloaded here [solution](https://github.com/bioinformaticsdotca/HTseq_2019/blob/master/IntegratedAssignment/integrative_assigment_commands.sh)



