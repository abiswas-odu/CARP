#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/M11
mkdir /scratch/mryan023/ISQuest/M11
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/M11_A_AGGCAGAA-TATCCTCT_L001.1.fastq /scratch/mryan023/FASTQ/M11_A_AGGCAGAA-TATCCTCT_L001.2.fastq /scratch/mryan023/ISQuest/ReadFiles/M11
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/M11/ /scratch/mryan023/ISQuest/ContigFiles/M11.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/M11.fasta M11 > /scratch/mryan023/ISQuest/M11/log.txt

