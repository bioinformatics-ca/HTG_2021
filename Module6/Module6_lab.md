---
layout: tutorial_page
permalink: /HTG_2021_module6_lab
title: HTG 2021 Module 6 Lab
header1: Workshop Pages for Students
header2: High-throughput Genomics Analysis Module 6 Lab
image: /site_images/CBW_High-throughput_icon.jpg
home: https://bioinformaticsdotca.github.io/HTG_2021
description: HTG 2021 Module 6 Lab
author: Jared Simpson
modified: September, 2021
---

# Genome Assembly for Short and Long Reads

by [Jared Simpson](https://simpsonlab.github.io)

## Introduction

In this lab we will perform de novo genome assembly of a bacterial genome using three different sequencing technologies. You will be guided through the genome assembly. At the end of the lab you will know:

1. How to run a short read assembler on Illumina data
2. How to run a long read assembler on Pacific Biosciences or Oxford Nanopore data
3. How to polish a long read assembly
4. How to assess the quality of an assembly

## Data Sets

In this lab we will use a bacterial data set to demonstrate genome assembly. This data set consists of sequencing reads for Escherichia coli. E. coli is a great data set for getting started on genome assembly as it is a small genome (4.6 Mbp) with relatively few repeats, and has a high quality reference. We have provided Illumina, PacBio HiFi and Oxford Nanopore reads for this genome. In this lab, you will make assemblies of the complete E. coli genome using `spades` and `flye`. 

## Data Preparation

First, lets create and move to a directory that we'll use to work on our assemblies:

```
mkdir -p ~/workspace/HTG/Module6
cd ~/workspace/HTG/Module6
```

For convenience, we'll make symbolic links to the data sets that we'll work with. Run this command from the terminal, which will find all of the sequencing data sets we provided (using the `ls` command) and symlink those files into your current working directory.


```
ls ~/CourseData/HTG_data/Module6/ecoli* | xargs -i ln -s {}
```

If you run `ls` you should now be able to see the three datasets. There are two files for the Illumina data as the two halves of the paired end reads are in separate files.

## E. coli Genome Assembly with Short Reads

Now we'll assemble the E. coli Illumina data using the [spades](http://bioinf.spbau.ru/spades) assembler. Parameterizing a short read assembly can be tricky and tuning the parameters (for example the size of the *k*-mer used) is often quite time consuming. Thankfully, spades will automatically select values for its parameters, making it particularly easy to use. You can start spades with this command (it will take a few minutes to run):

```
spades.py -o ecoli_illumina_spades/ -t 4 -1 ecoli_illumina_1.fastq -2 ecoli_illumina_2.fastq
```

After the assembly completes, let's move the results to a new directory that we'll use to keep track of all of our assemblies

```
mkdir -p assemblies
cp ecoli_illumina_spades/contigs.fasta assemblies/ecoli_illumina_spades.fasta
```

We can now start assessing the quality of our assembly. We typically measure the quality of an assembly using three factors:

- Contiguity: Long contigs are better than short contigs as long contigs give more information about the structure of the genome (for example, the order of genes)
- Completeness: Most of the genome should be assembled into contigs with few regions missing from the assembly
- Accuracy: The assembly should have few large-scale *misassemblies* and *consensus errors* (mismatches or insertions/deletions)

We'll use `abyss-fac` to calculate how contiguous our spades assembly is. Typically there will be a lot of short "leftover" contigs consisting of repetitive or low-complexity sequence, or reads with a very high error rate that could not be assembled. We don't want to include these in our statistics so we'll only use contigs that are at least 500bp in length (protip: piping tabular data into `column -t` will format the output so the columns nicely line up):

```
abyss-fac.pl -t 500 assemblies/ecoli_illumina_spades.fasta | column -t
```

The N50 statistic is the most commonly used measure of assembly contiguity. An N50 of x means that 50% of the assembly is represented in contigs x bp or longer. What is the N50 of the spades assembly? How many contigs were produced?

## E. coli Genome Assembly with Long Reads

Now, we'll use long sequencing reads to assemble the E. coli genome. Long sequencing reads are better at resolving repeats and typically give much more contiguous assemblies. Long reads need to use a different assembly strategy. For this tutorial, we'll use [flye](https://github.com/fenderglass/Flye) to assemble the PacBio HiFi dataset. Flye is a state-of-the-art assembler that is quite fast to run. Like Spades it is easy to configure as we only need to tell it the type of data we are using:

```
flye --out-dir ecoli_pacbio_flye --threads 4 --pacbio-hifi ecoli_pacbio.fastq
```

After the assembly completes, copy it to your aseemblies directory:

```
cp ecoli_pacbio_flye/assembly.fasta assemblies/ecoli_pacbio_flye.fasta
```

Our data set also includes an Oxford Nanopore data set. We can now assemble it using flye - the command is very similar, we just need to tell flye we are assembling nanopore data and change the input filename:

```
flye --out-dir ecoli_nanopore_flye --threads 4 --nano-raw ecoli_nanopore.fastq
```

And we'll copy the assembly:

```
cp ecoli_nanopore_flye/assembly.fasta assemblies/ecoli_nanopore_flye.fasta
```

## Assessing the Quality of your Assemblies using a Reference

The accuracy of the genome assembly is determined by how many misassemblies (large-scale rearrangements) and consensus errors (mismatches, insertions or deletions) the assembler makes. Calculating the accuracy of an assembly typically requires the use of a reference genome. We will use the [QUAST](http://quast.bioinf.spbau.ru/) software package to assess the accuracy of the assemblies.

Run QUAST on your three E. coli assemblies by running this command:

```
quast.py -R ~/CourseData/HTG_data/Module6/references/ecoli_k12.fasta assemblies/*.fasta
```

Using the web browser for your instance, open the QUAST PDF report (Module6/quast_results/latest/report.pdf) and try to determine which of the assemblies was a) the most complete b) the most contiguous and c) the most accurate.

## Assembly Polishing

The nanopore assembly can be improved by running a "polishing" step. There are many assembly polishing programs available for nanopore data (nanopolish, racon). To demonstrate polishing we will use a program called `medaka` that is particularly fast and easy to run. 

We're now going to use `medaka` to improve our assembly. Medaka uses a neural network which is trained to calculate a better consensus sequence for nanopore assemblies:

```
medaka_consensus -i ecoli_nanopore.fastq -d assemblies/ecoli_nanopore_flye.fasta -o ecoli_nanopore_medaka_polished -t 1
```

Now we can copy the medaka assembly to our output directory:

```
cp ecoli_nanopore_medaka_polished/consensus.fasta assemblies/ecoli_nanopore_flye_medaka.fasta
```

Now, re-run the QUAST step from above:

```
quast.py -R ~/CourseData/HTG_data/Module6/references/ecoli_k12.fasta assemblies/*.fasta
```

The report will be updated in Module6/quast_results/latest/report.html (all versions will also be stored in their own time-stamped directories in Module6/quast_results). Did the quality of your nanopore assembly improve?

## Bonus Exercise (time permitting)

Now, you can try to make a _hybrid_ assembly where reads from two technologies are assembled together. Try to use `spades` to make a hybrid assembly of the Illumina and nanopore datasets. We're not going to give you the command for this step, so you have a chance to run it on your own (hint: if you run `spades.py --help` it will print the list of options that the program takes). If you are stuck ask an instructor for help!

How does the hybrid assembly compare to the illumina-only or nanopore-only assemblies? 
