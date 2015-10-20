#!/bin/bash
#$ -N bwa 
#$ -M your.email@service.com
#$ -m beas
#$ -pe smp 2 
#$ -l h_vmem=80G
#$ -l h_rt=20:00:00
#$ -cwd
#$ -S /bin/bash

time ./01_scripts/01_bwa_w_fa.sh

