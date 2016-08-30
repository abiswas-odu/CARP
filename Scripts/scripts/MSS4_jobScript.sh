#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/MSS4
mkdir /scratch/mryan023/ISQuest/MSS4
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/MSS4_0708.1.fastq /scratch/mryan023/FASTQ/MSS4_0708.2.fastq /scratch/mryan023/ISQuest/ReadFiles/MSS4
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/MSS4/ /scratch/mryan023/ISQuest/ContigFiles/MSS4.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/MSS4.fasta MSS4 > /scratch/mryan023/ISQuest/MSS4/log.txt

