---
layout: tutorial_page
permalink: /HTSeq_2017_module6_lab_supplement
title: HT-Biology Lab 6
header1: Workshop Pages for Students
header2: Informatics on High-Throughput Sequencing Data
image: /site_images/CBW_High-throughput_icon.jpg
home: https://bioinformaticsdotca.github.io/htseq_2017
---

## Generating a preqc report for the E. coli data set

```
# First build an FM-index of the two E. coli Illumina data sets:
sga index -t 4 -a ropebwt ecoli.illumina.15x.fastq
sga index -t 4 -a ropebwt ecoli.illumina.50x.fastq

# Next, we can run `sga preqc` to run the calculations and generate the PDF report.
sga preqc -t 4 ecoli.illumina.15x.fastq > ecoli.illumina.15x.preqc
sga preqc -t 4 ecoli.illumina.50x.fastq > ecoli.illumina.50x.preqc
python sga-preqc-report.py ecoli.illumina.15x.preqc ecoli.illumina.50x.preqc
```

## Assembling the Oxford Nanopore Data

```
canu overlapper=mhap utgReAlign=true gnuplotTested=true -p ecoli-nanopore-canu -d ecoli-nanopore-auto genomeSize=4.6m -nanopore-raw ecoli.nanopore.50x.fastq
```


