#!/bin/bash
#SBATCH --job-name=nek_main
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --exclusive
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH --account=ebm5351_b_gpu
#SBATCH --partition=sla-prio

export CARDINAL_DIR=/storage/work/afc6440/david_cardinal_2/cardinal
export NEKRS_HOME=/storage/work/afc6440/david_cardinal_2/cardinal/install
export PYTHONPATH=/storage/work/afc6440/david_cardinal_2/cardinal/contrib/moose/python
export MOOSE_DIR=/storage/work/afc6440/david_cardinal_2/cardinal/contrib/moose
export LIBMESH_DIR=/storage/work/afc6440/david_cardinal_2/cardinal/contrib/moose/libmesh/installed
module purge
module load anaconda/2023.09
conda activate david-cardinal
module load openmpi/4.1.1-pmi2

module load gcc/9.1.0
export GCC_DIR=/storage/icds/RISE/sw8/gcc/gcc-9.1.0/bin/
export PATH=$GCC_DIR:$PATH 
export OMPI_CC=$GCC_DIR/gcc
export OMPI_CXX=$GCC_DIR/g++
export OMPI_FC=$GCC_DIR/gfortran

module load cuda/11.5.0
module load cmake/3.21.4

export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

mpirun -np 2 /storage/work/afc6440/david_cardinal_2/cardinal/cardinal-opt -i nek_main.i --nekrs-setup cylinder --nekrs-backend CUDA --nekrs-device-id 0
