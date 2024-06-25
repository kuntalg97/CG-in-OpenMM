#!/bin/bash

#SBATCH --job-name=total
#SBATCH --output=us-%j.out
#SBATCH --account=pi-gavoth
#SBATCH --partition=gavoth
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

#source activate /project2/gavoth/kuntalg/Softwares/openmm/openmm
conda activate /project2/gavoth/kuntalg/Softwares/PBNN/pbnnenv/

python CG_simulate.py
