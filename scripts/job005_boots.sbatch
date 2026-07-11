#!/bin/bash -l

#SBATCH --job-name=alf_boot
#SBATCH --time=0-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --array=1-100
#SBATCH --mail-user=ftsiboe@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

module purge
module load R/4.2.1-foss-2022a GDAL/3.5.0-foss-2022a PROJ/9.0.0-GCCcore-11.3.0
R --no-save -q < scripts/005_WarmingImpactsAlfalfa_boots.R
