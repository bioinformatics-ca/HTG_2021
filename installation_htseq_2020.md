---
layout: tutorial_page
permalink: /htseq_2020_installation
title: HTseq Installation Instructions
header1: Workshop Pages for Students
header2: Informatics for High-throughput Sequencing Data Analysis 2020 Installation Instructions
image: /site_images/CBW_High-throughput_icon.jpg
home: https://bioinformaticsdotca.github.io/htseq_2020
---

## Installation Instructions

1) Install latest version of R which can be downloaded from http://probability.ca/cran/.

1b) Download and install the most recent version of [R Studio desktop](http://www.rstudio.com/).  If prompted to install git, select yes.

2) Install the BioConductor core packages. If you have installed R version 3.5.0 or higher, open R and at the '>' prompt, paste the commands:
 
```
install.packages("BiocManager");
library(BiocManager);
BiocManager::install();
```

If you already have an older version of R installed (3.4.4 or lower), open R and at the '>' prompt, paste the commands:

```
source("http://bioconductor.org/biocLite.R");
biocLite();
```

If you are unsure which version you have installed, open R and at the '>' prompt, enter the command:

```
version
```

3) A robust text editor.   

* For Windows/PC - [notepad++](http://notepad-plus-plus.org/)  
* For Linux - [gEdit](http://projects.gnome.org/gedit/)  
* For Mac – [TextWrangler](http://www.barebones.com/products/textwrangler/download.html)

4) A file decompression tool.  

* For Windows/PC – [7zip](http://www.7-zip.org/).  
* For Linux – [gzip](http://www.gzip.org).   
* For Mac – already there.

5) A robust internet browser such as Firefox or Safari (Internet Explorer and Chrome are not recommended because of Java issues).

6) Java -The visualization program that we will be using (IGV) requires Java. Check if you have Java installed: https://www.java.com/verify/ and download Java if you do not have it installed (You need Java 11).

7) Integrative Genomics Viewer 2.4.x (IGV) - Once java is installed, go to https://software.broadinstitute.org/software/igv/download Click on the appropriate launch button that matches your computer's operating system.   

8) SSH client - Mac and Linux users already have a command line ssh program that can be run from the terminal. For Windows users, please download [PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html).  

9) SCP/SFTP client - We will be moving data from the servers to the student laptops for visualization. Mac and Linux users already have a command line scp and sftp program. For Windows users, please install [WinSCP](http://winscp.net/eng/download.php).

10) A PDF viewer (Adobe Acrobat or equivalent).
