#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/M13
mkdir /scratch/mryan023/ISQuest/M13
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/M13_small_0708.075x.1.fastq /scratch/mryan023/FASTQ/M13_small_0708.075x.2.fastq /scratch/mryan023/ISQuest/ReadFiles/M13
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/M13/ /scratch/mryan023/ISQuest/ContigFiles/M13.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/M13.fasta M13 > /scratch/mryan023/ISQuest/M13/log.txt

