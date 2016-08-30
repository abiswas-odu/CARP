#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/Rp72a
mkdir /scratch/mryan023/ISQuest/Rp72a
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/Rp72a_A_CGAGGCTG-CTCTCTAT_L001.1.fastq /scratch/mryan023/FASTQ/Rp72a_A_CGAGGCTG-CTCTCTAT_L001.2.fastq /scratch/mryan023/ISQuest/ReadFiles/Rp72a
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/Rp72a/ /scratch/mryan023/ISQuest/ContigFiles/Rp72a.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/Rp72a.fasta Rp72a > /scratch/mryan023/ISQuest/Rp72a/log.txt

