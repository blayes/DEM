#!/bin/bash
#$ -N sim_dem
#$ -q LT
#$ -l mf=32G
#$ -pe smp 2
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/dem/mixef/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-40
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R

R CMD BATCH --no-save --no-restore "--args 0 $SGE_TASK_ID" mixef_submit_dem.R sim_$SGE_TASK_ID.rout
