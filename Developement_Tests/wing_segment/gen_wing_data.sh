#!/bin/bash

##SBATCH --partition=compute
#SBATCH --partition=gpu
#SBATCH --qos=normal

##SBATCH --output=/scratch/msantana/%x-%N-%j.out # Output file
##SBATCH --error=/scratch/msantana/%x-%N-%j.err  # Error File

#SBATCH --time=3-00:00:00                 # run time limit (DD-HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --array=1
#SBATCH --ntasks=1		  # This for multiple workers for parallel distributed jobs
##SBATCH --cpus-per-task=56
#SBATCH --cpus-per-task=128                # This for multi threading. How many CPUS to use


module load matlab

#matlab -nodisplay -nosplash -nodesktop -nojvm -r "wing_section_self_convergence"
matlab -nodisplay -nosplash -nodesktop -nojvm -r "wing_section_pole_convergence"


