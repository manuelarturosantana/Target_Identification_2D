#!/bin/bash

#SBATCH --job-name=inspect_plane_data_360
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# wavefield-02-01 hangs on MATLAB batch jobs -- see plot_pole_eigenfunction.sh.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/inspect_plane_data_360_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/inspect_plane_data_360_%j.err

module load matlab

matlab -nodisplay -nosplash -nodesktop -batch "inspect_plane_data_360"

exit_code=$?
exit "$exit_code"
