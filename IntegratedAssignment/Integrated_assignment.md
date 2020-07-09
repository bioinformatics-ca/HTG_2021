---
layout: tutorial_page
permalink: /htseq_2020_IA
title: HTseq Integrated Assignment
header1: Workshop Pages for Students
header2: Informatics for High-throughput Sequencing Data Analysis 2020 Integrated Assignment
image: /site_images/CBW_High-throughput_icon.jpg
home: https://bioinformaticsdotca.github.io/htseq_2020
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

### Accessing a working node

When you log into the server, you are assigned to a "login" node (sometimes called a "head node"), which is shared by other users who are also logged in. As these nodes are a shared resouces, running computationally heavy workloads here can make the system unstable for everybody. In order to run your analysis in a stable environment without affecting other user you need to access a work node (sometimes called a "compute node"). Usually each job shoule be launched through the scheduler to run in a working environment, but our jobs in this workshop as are small and fast, so we can instead launch an interactive session on one of the work nodes by running:

```bash
salloc --mem 0 -n 8
```

The salloc command will assign us to a compute node and give us permission to use up to 8 cpus at a time. The interactive session will last for 1h, after which our session will end and we will be returned to the login node.


```bash
#set up environment variables
export WORK_DIR=$HOME/workspace/HTseq/Integrative_Assignment
export REF=$WORK_DIR/reference

rm -rf $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR
ln -fs ~/CourseData/HT_data/Module3/* .


# Load the software modules
module load \
  mugqic/java/openjdk-jdk1.8.0_72 \
  mugqic/bvatools/1.6 \
  mugqic/trimmomatic/0.36 \
  mugqic/samtools/1.9 \
  mugqic/bwa/0.7.17 \
  mugqic/GenomeAnalysisTK/4.1.0.0 \
  mugqic/R_Bioconductor/3.5.0_3.7
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

