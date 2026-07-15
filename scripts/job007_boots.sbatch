#!/bin/bash -l

#SBATCH --job-name=alf_boot
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --array=1-101
#SBATCH --mail-user=ftsiboe@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

# One array task = one bootstrap ("0000" + 100 draws => 101 tasks).
# Within a task the 50 GW specs (10 gwkit distance presets x 5 kernels) are run
# in parallel across --cpus-per-task cores (via future::multicore) and reduced to
# a per-county consensus in-task. Memory scales with the number of concurrent
# specs, not the number of cells; raise --mem/--cpus-per-task together.
#
# BEFORE submitting the array, run the one-off bandwidth cache (uncomment the
# "Bandwidth cache" function block in the boots script and source it once) so the
# array reuses fixed bandwidths instead of cross-validating every boot.

module purge
module load RStudio-Server/2023.09.1+494-foss-2023a-Java-11-R-4.3.2
R --no-save -q < scripts/007_WarmingImpactsAlfalfa_boots.R
