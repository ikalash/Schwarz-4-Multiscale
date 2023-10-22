#!/bin/bash
#SBATCH --job-name=schwarz
#SBATCH --time=360:00
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --mem=16GB
#SBATCH -o /scratch/users/jinnyc/test/optimize/%j.out

module load system
module load git
module load python/3.9.0
python3 $HOME/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/pinnschwarz/optimizer.py > /scratch/users/jinnyc/test/optimize/102223_run1.txt
