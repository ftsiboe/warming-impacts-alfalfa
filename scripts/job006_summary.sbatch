#!/bin/bash -l

#SBATCH --job-name=Alf_sum
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --array=1-45
#SBATCH --mail-user=ftsiboe@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

module purge
module load R/4.2.1-foss-2022a GDAL/3.5.0-foss-2022a PROJ/9.0.0-GCCcore-11.3.0
R --no-save -q < scripts/006_WarmingImpactsAlfalfa_summary.R
