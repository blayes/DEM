#!/bin/bash
#$ -N ml_iem
#$ -q INFORMATICS
#$ -pe 16cpn 32
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/dem/ml/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/
#$ -j y

tmphosts=`mktemp`
awk '{ for (i=0; i < $2; ++i) { print $1} }' $PE_HOSTFILE > $tmphosts

module load openmpi/intel-composer_xe_2015.3.187-1.6.5
eval "export $(mpirun env | grep OMPI_MCA_orte_precondition_transports)"
module load R/intel-composer_xe_2015.3.187_3.3.0

echo "Got $NSLOTS slots"
echo "jobid $JOB_ID"

mpirun -np 21 -machinefile $tmphosts R CMD BATCH --no-save --no-restore "--args 20 $SGE_TASK_ID" ml_iem_mpi.R iem_$SGE_TASK_ID.rout

exit 0
