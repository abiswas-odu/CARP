#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/1218R

mkdir /scratch/mryan023/ISQuest/1218R
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/1218R_A_TAAGGCGA-TAGATCGC_L001.1.fastq /scratch/mryan023/FASTQ/1218R_A_TAAGGCGA-TAGATCGC_L001.2.fastq /scratch/mryan023/ISQuest/ReadFiles/1218R
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/1218R/ /scratch/mryan023/ISQuest/ContigFiles/1218R.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/1218R.fasta 1218R > /scratch/mryan023/ISQuest/1218R/log.txt

