#!/bin/bash
#SBATCH --job-name=testSchwarz
#SBATCH --time=360:00
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --mem=16GB
module load system
module load git
module load python/3.9.0
python3 $HOME/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/pinnschwarz/opt_schwarz.py > /scratch/users/siqima/test_optimizer/run_schedule.txt
