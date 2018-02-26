#!/bin/bash
mkdir chim
echo ">>> aligning example reads"
# align RNA-seq reads to genome reference with bowtie2
# /project/cic3/Genome/hg_ncRNA/hg_All -- bowtie2 index of human, file.fastq -- your raw RNA-seq file, file.sam -- bowtie2 output file
# bowtie2 -p16 --very-sensitive --mm -D20 --score-min=C,-15,0 --phred33 -x /project/cic3/Genome/hg_ncRNA/hg_All -U file.fastq -S ./chim/file.sam
bowtie2 -p16 --very-sensitive --mm -D20 --score-min=C,-15,0 --phred33 -x /project/cic3/Genome/hg_ncRNA/hg_All -U /project/cic5/pso_GSE67785/PN/SRR1970400.fastq.bz2 -S ./chim/file.sam
# -p processing num?
#
# --very-sensitive Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
# -D puts an upper limit on the number of dynamic programming problems (i.e. seed extensions) that can “fail” in a row before Bowtie 2 stops searching. 
## Increasing -D makes Bowtie 2 slower, but increases the likelihood that it will report the correct alignment for a read that aligns many places.
# -R sets the maximum number of times Bowtie 2 will “re-seed” when attempting to align a read with repetitive seeds. 
## Increasing -R makes Bowtie 2 slower, but increases the likelihood that it will report the correct alignment for a read that aligns many places.
# -N Sets the number of mismatches to allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1. 
## Setting this higher makes alignment slower (often much slower) but increases sensitivity. 
### Default: 0.
# -L Sets the length of the seed substrings to align during multiseed alignment. 
## Smaller values make alignment slower but more sensitive. 
### Default: the --sensitive preset is used by default, which sets -L to 20 both in --end-to-end mode and in --local mode.
# -i Sets a function governing the interval between seed substrings to use during multiseed alignment. 
## For instance, if the read has 30 characters, and seed length is 10, and the seed interval is 6, the seeds extracted will be:
### Read:      TAGCTACGCTCTACGCTATCATGCATAAAC
### Seed 1 fw: TAGCTACGCT
### Seed 1 rc: AGCGTAGCTA
### Seed 2 fw:       CGCTCTACGC
### Seed 2 rc:       GCGTAGAGCG
### Seed 3 fw:             ACGCTATCAT
### Seed 3 rc:             ATGATAGCGT
### Seed 4 fw:                   TCATGCATAA
### Seed 4 rc:                   TTATGCATGA
### Since it’s best to use longer intervals for longer reads, 
### this parameter sets the interval as a function of the read length, 
### rather than a single one-size-fits-all number. 
#### For instance, specifying -i S,1,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length. 
#### See also: setting function options. If the function returns a result less than 1, it is rounded up to 1. 
#### Default: the --sensitive preset is used by default, which sets -i to S,1,1.15 in --end-to-end mode to -i S,1,0.75 in --local mode.
#
# -mm Use memory-mapped I/O to load the index, rather than typical file I/O. 
# Memory-mapping allows many concurrent bowtie processes on the same computer to share the same memory image of the index
# (i.e. you pay the memory overhead just once). 
# This facilitates memory-efficient parallelization of bowtie in situations where using -p is not possible or not preferable.
#
# -D20 
# 
# --score-min <func> Sets a function governing the minimum alignment score needed for an alignment to be considered “valid” (i.e. good enough to report).
# This is a function of read length. 
# For instance, specifying L,0,-0.6 sets the minimum-score function f to f(x) = 0 + -0.6 * x, where x is the read length. 
# See also: setting function options. The default in --end-to-end mode is L,-0.6,-0.6 and the default in --local mode is G,20,8.
#
# --score-min=C
#
# --phred33
## Input qualities are ASCII chars equal to the Phred quality plus 33. 
## This is also called the “Phred+33” encoding, which is used by the very latest Illumina pipelines.
# --phred64
## Input qualities are ASCII chars equal to the Phred quality plus 64. 
## This is also called the “Phred+64” encoding.
#
# -x <bt2-idx>
# The basename of the index for the reference genome. 
# The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. 
# bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.
#
# -U Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. 
# Reads may be a mix of different lengths. If - is specified, bowtie2 gets the reads from the “standard in” or “stdin” filehandle.
#
# -S <sam> File to write SAM alignments to. 
# By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console).
### bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i>} -S [<sam>]



echo ">>> get the unmapped"
# convert sam file to compressed bam file and sort it
samtools view -hbuS ./chim/file.sam | samtools sort - ./chim/file
# The -b, -C, -1, -u, -h, -H, and -c options change the output format from the default of headerless SAM, 
# and the -o and -U options set the output file name(s).
# 

# extract the unmapped reads from the sorted bam file
samtools view -hf 4 ./chim/file.bam | samtools view -Sb - > ./chim/unmapped_file.bam
rm ./chim/file.sam



echo ">>> split into anchors"
# get anchors (20nt each end of the read) from the unmapped reads and output a fasta(qfa) file
/project/cic3/xiaox/unmapped2anchors.py ./chim/unmapped_file.bam ./chim/anchors_file.qfa
gzip ./chim/anchors_file.qfa

echo ">>> run find_chimeric.py"
mkdir ./chim/file_chim
# map the anchors to the reference genome
bowtie2 --reorder --mm -D20 --score-min=C,-15,0 -q -x /project/cic3/Genome/GRCh38/GRCh38 -U ./chim/anchors_file.qfa.gz -S ./chim/file_candidate.sam
# find chimerich RNAs from the alignment results
python /project/cic3/xiaox/find_chimeric.py ./chim/file_candidate.sam ./chim/file_chim/sites.bed ./chim/file_chim/sites.reads -G /project/cic3/Genome/GRCh38/chromosome -p file_test -s ./chim/file_chim/sites.log -g /project/cic3/Genome/GRCh38/hg38/hg38_gene_sorted.bed.gz -d 0

