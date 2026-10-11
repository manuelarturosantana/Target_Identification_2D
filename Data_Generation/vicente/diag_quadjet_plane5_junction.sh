#!/bin/bash

#SBATCH --job-name=diag_quadjet_junction
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# wavefield-02-01 hangs on MATLAB batch jobs -- see plot_pole_eigenfunction.sh.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/diag_quadjet_junction_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/diag_quadjet_junction_%j.err

module load matlab

matlab -nodisplay -nosplash -nodesktop -batch "diag_quadjet_plane5_junction"

exit_code=$?
exit "$exit_code"
