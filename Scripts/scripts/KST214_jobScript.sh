#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/KST214
mkdir /scratch/mryan023/ISQuest/KST214
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/KST214_0708.050x.1.fastq /scratch/mryan023/FASTQ/KST214_0708.050x.2.fastq /scratch/mryan023/ISQuest/ReadFiles/KST214
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/KST214/ /scratch/mryan023/ISQuest/ContigFiles/KST214.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/KST214.fasta KST214 > /scratch/mryan023/ISQuest/KST214/log.txt

