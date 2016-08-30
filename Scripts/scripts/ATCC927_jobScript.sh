#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/ATCC927
mkdir /scratch/mryan023/ISQuest/ATCC927
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/ATCC927_2_ATGAGC_L001.050x.1.fastq /scratch/mryan023/FASTQ/ATCC927_2_ATGAGC_L001.050x.2.fastq /scratch/mryan023/ISQuest/ReadFiles/ATCC927
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/ATCC927/ /scratch/mryan023/ISQuest/ContigFiles/ATCC927.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/ATCC927.fasta ATCC927 > /scratch/mryan023/ISQuest/ATCC927/log.txt

