#!/bin/bash -l

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=blauch@cmu.edu
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0
#SBATCH --mem=2G

cd ~/git/familiarity_sims
conda activate mne

echo "$COMMAND"
$COMMAND
