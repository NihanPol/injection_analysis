#!/bin/sh

#PBS -l nodes=1:ppn=1,pvmem=50gb
##PBS -m abe
#PBS -M noemail@hpc.wvu.edu
##PBS -q stmcwilliams_lp
#PBS -q debug

#PBS -N injecting_real_data

module load compilers/gcc/8.2.0 compilers/python/2.7.13 libraries/suitesparse/5.3_gcc82 compilers/cuda/7.5
module load libraries/openblas/0.3.2_gcc82 mpi/openmpi/3.1.2_gcc82

export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd /users/nspol/stochastic_11yr_analysis/notebooks/new_injection_analysis/

python final_inject_gwb.py -amps_path injected_amps.npy -outdir ./test/ -seed 123 456
