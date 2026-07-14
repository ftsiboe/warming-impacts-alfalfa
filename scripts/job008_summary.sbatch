#!/bin/bash -l

#SBATCH --job-name=Alf_sum
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --mail-user=ftsiboe@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

# Single sequential job: the rewritten summary reads all per-cell consensus
# files and writes the summary_*.rds outputs. No SLURM array (it no longer
# partitions by GW spec).

module purge
module load RStudio-Server/2023.09.1+494-foss-2023a-Java-11-R-4.3.2
R --no-save -q < scripts/008_WarmingImpactsAlfalfa_summary.R
