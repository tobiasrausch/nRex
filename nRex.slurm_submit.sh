#!/bin/bash
#SBATCH -A korbel                   # group to which you belong
#SBATCH -p 1month                   # partition (queue)
#SBATCH -N 1                        # number of nodes
#SBATCH -n 4                        # number of cores
#SBATCH --mem 8000M                 # memory pool for all cores
#SBATCH -t 16-8:00                  # time (D-HH:MM)
#SBATCH -o nRex.%N.%j.out           # STDOUT
#SBATCH -e nRex.%N.%j.err           # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=rausch@embl.de  # send-to address

# Fetch ATAC-Seq script
NREXSCRIPT=${1}
shift

# Run analysis pipeline
${NREXSCRIPT} $@
