#!/bin/bash

#SBATCH --job-name=total
#SBATCH --output=us-%j.out
#SBATCH --account=pi-gavoth
#SBATCH --partition=gavoth
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load python

source activate /project2/gavoth/kuntalg/Softwares/openmm/openmm

python simulatewater.py
