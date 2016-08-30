#!/bin/tcsh

module load openmpi/open64/64/1.6.5
module load blast/2.2.29
module load bioperl/1.6.923
unlimit

rm -r /scratch/mryan023/ISQuest/453-022
mkdir /scratch/mryan023/ISQuest/453-022
/scratch/abisw001/ISQuest/fastQ2A.pl /scratch/mryan023/FASTQ/453-022_A_CAGAGAGG-CTAAGCCT_L001.100x.1.fastq /scratch/mryan023/FASTQ/453-022_A_CAGAGAGG-CTAAGCCT_L001.100x.2.fastq /scratch/mryan023/ISQuest/ReadFiles/453-022
/scratch/mryan023/ISQuest/isquest-code/isq /scratch/mryan023/ISQuest/MMAR32Param.conf /scratch/mryan023/ISQuest/453-022/ /scratch/mryan023/ISQuest/ContigFiles/453-022.cat.fasta /scratch/mryan023/ISQuest/ReadFiles/453-022.fasta 453-022 > /scratch/mryan023/ISQuest/453-022/log.txt

