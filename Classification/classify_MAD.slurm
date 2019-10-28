#!/bin/bash

outer_cv=20 # number of outer cross-validation iterations to perform. larger number -> less variance in results, greater power due to less variance in null distribution
overwrite='--overwrite' # set to '' to not overwrite existing results
nd='--compute-nulldist' # set to '' to avoid null dist computation (much faster)
declare -a mad_types=('MAD')
declare -a dtypes=('hit' 'FA' 'miss')
declare -a alts=('' '--alt' '--distractor-diff')
use_slurm=1 # set to 0 to run locally (will take some time due to null distribution calculation)

for mad_type in ${mad_types[@]}; do
  for dtype in ${dtypes[@]}; do
    for alt in "${alts[@]}"; do
      full_log=${log_name}-${dtype}${alt}${irf}-${mad_type}
      COMMAND="python classify_MAD.py --dtype $dtype $alt_tag $nd --outer-cv $outer_cv --measure $mad_type $overwrite"
      if [ ${use_slurm} -eq 1 ]; then
        sbatch --export="COMMAND=$COMMAND" --time 48:00:00 --output=log/slurm-%j.log run_slurm.sbatch
      else; then
        $COMMAND
      fi
    done
  done
done
