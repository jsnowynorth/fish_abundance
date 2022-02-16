#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH --job-name=spatial_fish_model  # name for the job
#SBATCH --cpus-per-task=1              # number of cores
#SBATCH --mem=100G                     # total memory
#SBATCH --time 1-20:00                 # time limit in the form days-hours:minutes
#SBATCH --mail-user=joshuanorth@mail.missouri.edu
#SBATCH --mail-type=FAIL,END           # email types
#SBATCH -A stsn                        # denotes part of space time group
#SBATCH --partition Lewis
#--------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

## Module Commands
# 'use module avail r/' to find the latest version
module load r
module list

## Run the R script
SCRIPT='lewis_spatial_no_land_script.R'
srun Rscript ${SCRIPT}

echo "### Ending at: $(date) ###"
