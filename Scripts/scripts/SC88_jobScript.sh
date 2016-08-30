#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/SC88
mkdir /scratch/mryan023/ISQuest/SC88
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/SC88_CTATAC_L001.075x.1.fastq /scratch/mryan023/FASTQ/SC88_CTATAC_L001.075x.2.fastq /scratch/mryan023/ISQuest/ReadFiles/SC88
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/SC88/ /scratch/mryan023/ISQuest/ContigFiles/SC88.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/SC88.fasta SC88 > /scratch/mryan023/ISQuest/SC88/log.txt

