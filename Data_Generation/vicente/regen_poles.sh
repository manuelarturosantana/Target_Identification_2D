#!/bin/bash

#SBATCH --job-name=regen_poles
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G

# job 13113 OOM'd (7 oom_kill events, ~6 min in) with 128 'Processes'
# workers each holding their own copy of the N=3520 matrix machinery --
# 32 workers is within the range that ran cleanly in diag_pole.sh (16
# workers, 32G, ~60 dense N=3520 solves with headroom to spare).

# wavefield-02-01 hangs on MATLAB batch jobs -- see plot_pole_eigenfunction.sh.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/regen_poles_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/regen_poles_%j.err

module load matlab

echo "Starting regen_poles config ${CONFIG_ID:-<default: 6, twinjet_plane1>}"

matlab -nodisplay -nosplash -nodesktop -batch "regen_poles"

exit_code=$?
exit "$exit_code"
