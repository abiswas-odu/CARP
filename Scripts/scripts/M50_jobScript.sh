#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/M50
mkdir /scratch/mryan023/ISQuest/M50
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/M50_0708.1.fastq /scratch/mryan023/FASTQ/M50_0708.2.fastq /scratch/mryan023/ISQuest/ReadFiles/M50
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/M50/ /scratch/mryan023/ISQuest/ContigFiles/M50.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/M50.fasta M50 > /scratch/mryan023/ISQuest/M50/log.txt

