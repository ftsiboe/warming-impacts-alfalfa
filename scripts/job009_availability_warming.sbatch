#!/bin/bash -l

#SBATCH --job-name=alf_avail
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60G
#SBATCH --mail-user=ftsiboe@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

# Single sequential job. Warming-scenario alfalfa availability + cattle shifts,
# computed OUTSIDE the bootstrap. Per analysis cell, the neighbourhood lag is
# applied ONCE per distance metric to EVERY boot x scenario production vector at
# the same time (estimate_gwlag_kernels reuses one distance matrix across the 5
# kernels and lags all columns in a single BLAS pass), so the GW work collapses
# from 50 specs x 101 boots to 10 distance builds per cell. Memory is dominated
# by the obs x obs weight matrix plus the wide value matrix; raise --mem if the
# county set or boot count grows.

module purge
module load RStudio-Server/2023.09.1+494-foss-2023a-Java-11-R-4.3.2
R --no-save -q < scripts/009_WarmingImpactsAlfalfa_availability_warming.R
