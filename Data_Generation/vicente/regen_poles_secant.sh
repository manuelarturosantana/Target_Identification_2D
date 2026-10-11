#!/bin/bash

#SBATCH --job-name=regen_poles_secant
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# wavefield-02-01 hangs on MATLAB batch jobs -- see plot_pole_eigenfunction.sh.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/regen_poles_secant_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/regen_poles_secant_%j.err

module load matlab

matlab -nodisplay -nosplash -nodesktop -batch "regen_poles_secant"

exit_code=$?
exit "$exit_code"
