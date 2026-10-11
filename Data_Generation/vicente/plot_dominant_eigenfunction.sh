#!/bin/bash

#SBATCH --job-name=dominant_eigfun
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G

# wavefield-02-01 hangs on MATLAB batch jobs -- see plot_pole_eigenfunction.sh.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/dominant_eigfun_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/dominant_eigfun_%j.err

module load matlab

echo "Config ${CONFIG_ID:-<default: 6, twinjet_plane1>}, ng=${NG:-40}"

matlab -nodisplay -nosplash -nodesktop -batch "plot_dominant_eigenfunction"

exit_code=$?
exit "$exit_code"
