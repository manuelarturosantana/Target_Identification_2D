#!/bin/bash

#SBATCH --job-name=replot_glrtvote_logscale
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# wavefield-02-01 hangs on MATLAB batch jobs -- see plot_pole_eigenfunction.sh.
#SBATCH --exclude=wavefield-02-01

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/replot_glrtvote_logscale_%j.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/replot_glrtvote_logscale_%j.err

module load matlab

matlab -nodisplay -nosplash -nodesktop -batch "replot_glrtvote_logscale"

exit_code=$?
exit "$exit_code"
