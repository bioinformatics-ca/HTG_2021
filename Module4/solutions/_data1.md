A bai file isn't an indexed form of a bam - it's a companion to your bam that contains the index.

A bam file is a binary blob that stores all of your aligned sequence data. The index allows to acces specific region of the bam file without having to decompress the data and read it from start. 


Working with indexed bam file allows to speed-up the data acces. Moreover many tools require index bam files.  
