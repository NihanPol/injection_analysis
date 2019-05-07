#!/bin/sh

#PBS -l nodes=1:ppn=2,pvmem=15gb
##PBS -m abe
#PBS -M noemail@hpc.wvu.edu
##PBS -q stmcwilliams_lp
#PBS -q comm_mmem_week

#PBS -t 20-29

#PBS -N injection_test_${PBS_ARRAYID}

module load compilers/python/2.7.13 environment/gcc-mpi mpi/openmpi/2.0.1 compilers/cuda/7.5 
module load compilers/gcc/7.3.0

cd /users/nspol/stochastic_11yr_analysis/notebooks/injection_analysis/

##The first input to the python script is realization number [0,9] while the second is the index
##of the amplitude array [0, 29].

python 11yr_injection_BE_3A.py 2 ${PBS_ARRAYID}
