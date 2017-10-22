#!/bin/bash
#$ -N lmer_ml
#$ -q INFORMATICS
#$ -pe 16cpn 16
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/dem/ml/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load openmpi/intel-composer_xe_2015.3.187-1.6.5
eval "export $(mpirun env | grep OMPI_MCA_orte_precondition_transports)"
module load R/intel-composer_xe_2015.3.187_3.3.0

R CMD BATCH --no-save --no-restore "--args 1 $SGE_TASK_ID" ml_submit_dem.R lmer/lmer_$SGE_TASK_ID.rout
