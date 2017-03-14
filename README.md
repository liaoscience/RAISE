######################RAISE is used to circRNA ReAlign Inner Splice and Expression (RAISE) for circRNA backsplice sites identification, quantification of abundance and prediction of their internal structure using RNA-seq data. 

######################


RAISE
==============

Contains the following files:
- pipeline_prepare.txt
- pipeline_4tools_align.bash
- pipeline_main_abundance.bash
- acfs_example.txt
- put_together_ws.pl

four circRNA detection tools
- find_circ
- mapsplice (2.1.8)
- acfs
- circRNA_finder


These scripts have been tested on various Linux distributions. Before they can be run, make sure that the following prerequisites are installed:
 - perl
 - awk
 - STAR (versions 2.3.1)
 - samtools (1.3)
 - bwa (0.7.3)
 - hisat2 (2.0.4)
 - bowtie2 (2.2.5)
 - bedtools (2.24.0)
 - stringtie (1.3.0)
 - htseq-count (0.6.0)


To run the scripts to identify circular RNAs, first run hisat2 and STAR, once for each data set:

sh pipeline_4tools_align.bash

Next, run the post processing scripts:

sh pipeline_main_abundance.bash

For each library the following output files or folds are produced:

folds:
a) <lib_name>_fc: find_circ output
b) <lib_name>_acfs: acfs output
c) <lib_name>_mapsplice: mapsplice output
d) circRNA_validation: circRNA realignment details
files:
a) <lib name>_filteredJunctions.bed: circRNA_finder output
b) <lib name>_total.txt: circRNA abundance
b) <lib name>_total.txt: circRNA abundance detail


any suggestion and question are welcome. email: liao_science@163.com

