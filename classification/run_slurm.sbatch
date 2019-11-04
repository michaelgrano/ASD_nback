#!/bin/bash -l

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=blauch@cmu.edu
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0
#SBATCH --mem=2G
#SBATCH --exclude=mind-0-13,mind-0-14

cd ~/git/nback_ASD/Classification
conda activate nback

echo "$COMMAND"
$COMMAND
