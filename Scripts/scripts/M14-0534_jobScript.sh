#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/14-0534
mkdir /scratch/mryan023/ISQuest/14-0534
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/14-0534_S3_L001.100x.1.fastq /scratch/mryan023/FASTQ/14-0534_S3_L001.100x.2.fastq /scratch/mryan023/ISQuest/ReadFiles/14-0534
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/14-0534/ /scratch/mryan023/ISQuest/ContigFiles/14-0534.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/14-0534.fasta 14-0534 > /scratch/mryan023/ISQuest/14-0534/log.txt

