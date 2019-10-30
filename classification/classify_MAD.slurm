#!/bin/bash

overwrite='--overwrite' # set to '' to not overwrite existing results
nd='' # set to '--no-nulldist' to avoid null dist computation (much faster)
declare -a dtypes=('hit' 'FA' 'miss')
declare -a dconds=('' '--distractors' '--distractor-diff')
use_slurm=0 # set to 0 to run locally (will take some time with null distribution calculation)

ii=0
for dtype in ${dtypes[@]}; do
  for distractors in "${dconds[@]}"; do
    let ii+=1
    COMMAND="python classify_MAD.py --dtype $dtype $distractors $nd $overwrite"
    if [ ${use_slurm} -eq 1 ]; then
      sbatch --export="COMMAND=$COMMAND" --job-name ASD-${ii} --time 48:00:00 --output=log/slurm-%j.log run_slurm.sbatch
    else
      $COMMAND
    fi
  done
done
