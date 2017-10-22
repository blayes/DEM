#!/bin/bash
#$ -N lmer_dem
#$ -q INFORMATICS
#$ -pe smp 4
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/dem/mixef/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-80
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R

R CMD BATCH --no-save --no-restore "--args 1 $SGE_TASK_ID" mixef_submit_dem.R lmer_$SGE_TASK_ID.rout
