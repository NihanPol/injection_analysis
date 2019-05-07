#!/bin/sh

#PBS -l nodes=1:ppn=2,pvmem=15gb
##PBS -m abe
#PBS -M noemail@hpc.wvu.edu
#PBS -q stmcwilliams_lp
##PBS -q comm_mmem_week

#PBS -t 0-29

#PBS -N OS_real_injection_2A_${PBS_ARRAYID}

module load compilers/gcc/8.2.0 compilers/python/2.7.13 libraries/suitesparse/5.3_gcc82 compilers/cuda/7.5
module load libraries/openblas/0.3.2_gcc82 mpi/openmpi/3.1.2_gcc82

#module load libraries/suitesparse/5.3_gcc82 libraries/openblas/0.3.2_gcc82 mpi/openmpi/3.1.2_gcc82

export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1 

cd /users/nspol/stochastic_11yr_analysis/notebooks/injection_analysis/

##The first input to the python script is realization number [0,9] while the second is the index
##of the amplitude array [0, 29].

python 11yr_OS_2A_noBE.py 0 ${PBS_ARRAYID}
