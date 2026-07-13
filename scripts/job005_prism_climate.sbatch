#!/bin/bash -l

#SBATCH --job-name=prism_climate1
#SBATCH --time=0-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=28G
#SBATCH --array=1-1000
#SBATCH --mail-user=ftsiboe@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

module purge
module load RStudio-Server/2023.09.1+494-foss-2023a-Java-11-R-4.3.2
R --no-save -q < scripts/005_WarmingImpactsAlfalfa_data_prism_climate.R

