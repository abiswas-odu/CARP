#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/ATCC29254
mkdir /scratch/mryan023/ISQuest/ATCC29254
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/ATCC29254_0708.1.fastq /scratch/mryan023/FASTQ/ATCC29254_0708.2.fastq /scratch/mryan023/ISQuest/ReadFiles/ATCC29254
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/ATCC29254/ /scratch/mryan023/ISQuest/ContigFiles/ATCC29254.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/ATCC29254.fasta ATCC29254 > /scratch/mryan023/ISQuest/ATCC29254/log.txt

