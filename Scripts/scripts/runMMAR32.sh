#!/bin/tcsh

qsub -pe openmpi 1 ./M1218R_jobScript.sh
qsub -pe openmpi 1 ./M12-8618_jobScript.sh
qsub -pe openmpi 1 ./M13-3845_jobScript.sh
qsub -pe openmpi 1 ./M14-0534_jobScript.sh
qsub -pe openmpi 1 ./M324-958_jobScript.sh
qsub -pe openmpi 1 ./M453-022_jobScript.sh
qsub -pe openmpi 1 ./ATCC11564_jobScript.sh
qsub -pe openmpi 1 ./ATCC29254_jobScript.sh
qsub -pe openmpi 1 ./ATCC927_jobScript.sh
qsub -pe openmpi 1 ./C35_jobScript.sh
qsub -pe openmpi 1 ./C7_jobScript.sh
qsub -pe openmpi 1 ./Davis_jobScript.sh
qsub -pe openmpi 1 ./KST214_jobScript.sh
qsub -pe openmpi 1 ./KST266_jobScript.sh
qsub -pe openmpi 1 ./KST417_jobScript.sh
qsub -pe openmpi 1 ./KST687_jobScript.sh
qsub -pe openmpi 1 ./KST94_jobScript.sh
qsub -pe openmpi 1 ./M11_jobScript.sh
qsub -pe openmpi 1 ./M13_jobScript.sh
qsub -pe openmpi 1 ./M30_jobScript.sh
qsub -pe openmpi 1 ./M50_jobScript.sh
qsub -pe openmpi 1 ./MSS4_jobScript.sh
qsub -pe openmpi 1 ./R106_jobScript.sh
qsub -pe openmpi 1 ./R171_jobScript.sh
qsub -pe openmpi 1 ./Rp72a_jobScript.sh
qsub -pe openmpi 1 ./SC86_jobScript.sh
qsub -pe openmpi 1 ./SC88_jobScript.sh

