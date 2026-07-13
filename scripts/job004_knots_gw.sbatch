#!/bin/bash -l

#SBATCH --job-name=alf_knots_gw
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=80G
#SBATCH --mail-user=ftsiboe@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

# One-off, full-sample county degree-day knot search (NOT bootstrapped, NOT an
# array). The 50 GW specs run in parallel across --cpus-per-task cores; each spec
# sweeps the +/- knot_band candidate pairs and picks each county's best valid
# knot, then knots are consensused across specs -> output/optimal_knots_gw.rds.
# Run after 003 (national knot) and before 005 stage-2 (dday trim).

module purge
module load RStudio-Server/2023.09.1+494-foss-2023a-Java-11-R-4.3.2
R --no-save -q < scripts/004_WarmingImpactsAlfalfa_knots_gw.R
