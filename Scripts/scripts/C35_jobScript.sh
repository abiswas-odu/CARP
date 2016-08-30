#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/C35
mkdir /scratch/mryan023/ISQuest/C35
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/C35_0708.1.fastq /scratch/mryan023/FASTQ/C35_0708.2.fastq /scratch/mryan023/ISQuest/ReadFiles/C35
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/C35/ /scratch/mryan023/ISQuest/ContigFiles/C35.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/C35.fasta C35 > /scratch/mryan023/ISQuest/C35/log.txt

