#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/R106
mkdir /scratch/mryan023/ISQuest/R106
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/R106_A_CGTACTAG-CTCTCTAT_L001.1.fastq /scratch/mryan023/FASTQ/R106_A_CGTACTAG-CTCTCTAT_L001.2.fastq /scratch/mryan023/ISQuest/ReadFiles/R106
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/R106/ /scratch/mryan023/ISQuest/ContigFiles/R106.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/R106.fasta R106 > /scratch/mryan023/ISQuest/R106/log.txt

