#!/bin/bash
#SBATCH -A korbel                   # group to which you belong
#SBATCH -p 1month                   # partition (queue)
#SBATCH -N 1                        # number of nodes
#SBATCH -n 4                        # number of cores
#SBATCH --mem 16000M                # memory pool for all cores
#SBATCH -t 10-2:00                  # time (D-HH:MM)
#SBATCH -o nRex.%N.%j.out           # STDOUT
#SBATCH -e nRex.%N.%j.err           # STDERR
#SBATCH --mail-type=FAIL            # notifications for job done & fail
#SBATCH --mail-user=rausch@embl.de  # send-to address

# Do we have EasyBuild and module
module -v > /dev/null 2>&1 || { echo >&2 "EasyBuild modules are required. Aborting."; exit 1; }

# Load required modules
module load Java
module load SAMtools
module load BCFtools
module load bwa
module load BEDTools
module load R-bundle-Bioconductor
module load expat
module load Boost
module load HTSlib
module load alfred
module load Perl
module load Python/2.7.12-foss-2016b

# Fetch ATAC-Seq script
NREXSCRIPT=${1}
shift

# Run analysis pipeline
${NREXSCRIPT} $@
