#!/bin/sh

#PBS -l nodes=1:ppn=2,pvmem=15gb
##PBS -m abe
#PBS -M noemail@hpc.wvu.edu
#PBS -q stmcwilliams_lp
##PBS -q comm_mmem_week

#PBS -t 0-14

#PBS -N OS_ipta_2A_${PBS_ARRAYID}

module load compilers/gcc/8.2.0 compilers/python/2.7.13 libraries/suitesparse/5.3_gcc82 compilers/cuda/7.5
module load libraries/openblas/0.3.2_gcc82 mpi/openmpi/3.1.2_gcc82

#module load libraries/suitesparse/5.3_gcc82 libraries/openblas/0.3.2_gcc82 mpi/openmpi/3.1.2_gcc82

export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1 

cd /users/nspol/ipta_cut_dr2_injection_analysis/injection_analysis/

#Stuff to avoid IO bomb
rsync -av /users/nspol/.local/tempo2/share/tempo2 /dev/shm/
export TEMPO2=/dev/shm/tempo2

##The first input to the python script is realization number [0,9] while the second is the index
##of the amplitude array [0, 29].

python 11yr_OS_2A_noBE.py -realiz 0 -amp_index ${PBS_ARRAYID} -chainpath /scratch/nspol/ipta_cut_results/2a_wo_BE/ -outdir /scratch/nspol/my_ipta_cut_results/2a_wo_BE_OS/ --psrlist /scratch/nspol/IPTA/partim_cut_IPTA/psrlist_cut_IPTA.txt --timpath /scratch/nspol/IPTA/ --parpath /scratch/nspol/IPTA/partim_cut_IPTA/ --noisepath /scratch/nspol/IPTA/partim_cut_IPTA/ --amps_path ipta_amps.npy
