#!/bin/bash

#SBATCH --job-name=constellation
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128

# CHANGE IF CLONED: folder where the MATLAB script lives.
#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente

# CHANGE IF CLONED: folder where SLURM logs should be written.
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/constellations_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/constellations_%j.err

module load matlab

# CHANGE IF CLONED: MATLAB script name, without ".m".
matlab -nodisplay -nosplash -nodesktop -nojvm \
    -r "data_gen_constellation; exit"