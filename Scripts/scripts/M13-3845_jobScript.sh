#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/13-3845
mkdir /scratch/mryan023/ISQuest/13-3845
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/13-3845_S10_L001.075x.1.fastq /scratch/mryan023/FASTQ/13-3845_S10_L001.075x.2.fastq /scratch/mryan023/ISQuest/ReadFiles/13-3845
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/13-3845/ /scratch/mryan023/ISQuest/ContigFiles/13-3845.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/13-3845.fasta 13-3845 > /scratch/mryan023/ISQuest/13-3845/log.txt

