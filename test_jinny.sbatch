#!/bin/bash
#SBATCH --job-name=schwarz
#SBATCH --time=360:00
#SBATCH -p normal
#SBATCH --cpus-per-task=10
#SBATCH --mem=16GB
#SBATCH -o /scratch/users/jinnyc/test/optimize/%j.out

module load system
module load git
module load python/3.9.0
python3 $HOME/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/pinnschwarz/opt_PINN.py > /scratch/users/jinnyc/test/optimize/110623-3.txt
