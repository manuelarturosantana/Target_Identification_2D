#!/bin/bash

#SBATCH --job-name=late_time_nacelle_field
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=150G

# One array task per plane in the classification experiment (config ids
# 1-5 quadjet, 6-10 twinjet -- see the configs table in
# late_time_nacelle_field.m). Capped at 3 concurrent tasks so this batch
# doesn't monopolize the partition overnight; matches the %3 convention
# already used by deformation.sh/constellation.sh in this folder.
##SBATCH --array=1-10%3

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/late_time_nacelle_field_%A_%a.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/late_time_nacelle_field_%A_%a.err

module load matlab

status_dir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/status/late_time_nacelle_field
mkdir -p "$status_dir"

echo "Starting late_time_nacelle_field config ${SLURM_ARRAY_TASK_ID:-<none>}"

# No -nojvm/-singleCompThread: parpool needs Java, and this script
# renders + exports a figure (exportgraphics), unlike the other data_gen
# scripts in this folder.
matlab -nodisplay -nosplash -nodesktop -batch "late_time_nacelle_field"

exit_code=$?

if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    if [ "$exit_code" -eq 0 ]; then
        echo "DONE config ${SLURM_ARRAY_TASK_ID}" \
            > "${status_dir}/config_${SLURM_ARRAY_TASK_ID}.done"
    else
        echo "FAILED config ${SLURM_ARRAY_TASK_ID} with exit code ${exit_code}" \
            > "${status_dir}/config_${SLURM_ARRAY_TASK_ID}.failed"
    fi
fi

exit "$exit_code"
