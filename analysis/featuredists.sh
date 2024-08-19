#!/bin/bash
#SBATCH --account=def-ijc               # account associated with job
#SBATCH --nodes=1                       # number of node MUST be 1
#SBATCH --cpus-per-task=44              # number of processes
#SBATCH --mem=0                         # memory usage
#SBATCH --time=0-03:00                  # time (DD-HH:MM)
#SBATCH --mail-user=mondrus@ualberta.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL                 # Send an email in all cases (job started, job ended, job aborted)

module load gcc/9.3.0 r/4.0.2
export R_LIBS=~/local/R_libs/
R CMD BATCH --no-save --no-restore featuredists.R