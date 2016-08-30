#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/R171
mkdir /scratch/mryan023/ISQuest/R171
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/R171_A_GTAGAGGA-AGAGTAGA_L001.1.fastq /scratch/mryan023/FASTQ/R171_A_GTAGAGGA-AGAGTAGA_L001.2.fastq /scratch/mryan023/ISQuest/ReadFiles/R171
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/R171/ /scratch/mryan023/ISQuest/ContigFiles/R171.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/R171.fasta R171 > /scratch/mryan023/ISQuest/R171/log.txt

