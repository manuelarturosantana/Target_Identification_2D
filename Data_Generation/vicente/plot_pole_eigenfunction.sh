#!/bin/bash

#SBATCH --job-name=pole_eigenfunction
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=150G

# wavefield-02-01 hangs on MATLAB batch jobs (confirmed on job 13061:
# 2+ hours RUNNING with zero script output and near-idle CPULoad despite
# 128 CPUs allocated -- never got past MATLAB/parpool startup). Fine for
# non-MATLAB jobs, excluded here because this script is MATLAB.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/pole_eigenfunction_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/pole_eigenfunction_%j.err

module load matlab

echo "Starting plot_pole_eigenfunction config ${CONFIG_ID:-<default: 6, twinjet_plane1>}"

matlab -nodisplay -nosplash -nodesktop -batch "plot_pole_eigenfunction"

exit_code=$?
exit "$exit_code"
