#!/bin/bash -l

#SBATCH --job-name=alf_knots_gw
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --array=1-10
#SBATCH --mail-user=ftsiboe@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

# One-off, full-sample county degree-day knot search (NOT bootstrapped). The
# array is BY DISTANCE METRIC (10 tasks); each task fits all 5 kernels for its
# metric via estimate_gwr_kernels(), which builds the obs->target distance matrix
# ONCE per (metric, pair) and reuses it across kernels, and selects the bandwidth
# ONCE per metric (reused across pairs). Each task writes one file per spec
# (metric x kernel) to fastscratch/knots_gw/.
#
# Consensus across the 50 specs -> output/optimal_knots_gw.rds is built only once
# ALL 50 spec files exist, so after the array completes run the script ONCE more
# (no array) to aggregate:  R --no-save -q < scripts/004_WarmingImpactsAlfalfa_knots_gw.R
# Run after 003 (national knot) and before 005 stage-2 (dday trim).

module purge
module load RStudio-Server/2023.09.1+494-foss-2023a-Java-11-R-4.3.2
R --no-save -q < scripts/004_WarmingImpactsAlfalfa_knots_gw.R
