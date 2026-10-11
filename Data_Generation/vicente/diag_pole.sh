#!/bin/bash

#SBATCH --job-name=diag_pole
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# wavefield-02-01 hangs on MATLAB batch jobs -- see plot_pole_eigenfunction.sh.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/diag_pole_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/diag_pole_%j.err

module load matlab

matlab -nodisplay -nosplash -nodesktop -batch "run('/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/diag_pole.m')"

exit_code=$?
exit "$exit_code"
