#!/bin/bash

#SBATCH --job-name=dominant_eigfun_glrtvote
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=300G

#SBATCH --array=1-10

# wavefield-02-01 hangs on MATLAB batch jobs -- see plot_pole_eigenfunction.sh.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/dominant_eigfun_glrtvote_%A_%a.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/dominant_eigfun_glrtvote_%A_%a.err

module load matlab

export NG=1024
export SAVE_DIR=/scratch/vhojas/data_planes_2d/glrtvote_figures
mkdir -p "$SAVE_DIR"

echo "Config ${SLURM_ARRAY_TASK_ID}, ng=${NG}, save_dir=${SAVE_DIR}"

matlab -nodisplay -nosplash -nodesktop -batch "plot_dominant_eigenfunction_glrtvote"

exit_code=$?
exit "$exit_code"
