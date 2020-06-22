

log_name=aut
binspace=10
outer_cv=20
irf=''
overwrite='--overwrite'
declare -a mad_types=('MAD') # 'MaxAD' 'MeanAD' 'MADs')
declare -a dtypes=('hit' 'FA' 'miss')
declare -a alts=('' '--alt')
declare -a alts=('--distractor-diff')

for mad_type in ${mad_types[@]}; do
  for dtype in ${dtypes[@]}; do
    for alt in "${alts[@]}"; do
      full_log=${log_name}-${dtype}${alt}${irf}-${mad_type}
      qsub classify_MAD.sh \
      -d /home/nblauch/git/autism \
      -N $full_log \
      -v binspace=${binspace},outer_cv=${outer_cv},dtype=${dtype},mad_type=${mad_type},overwrite="${overwrite}",alt_tag=${alt} \
      -l nodes=1:ppn=1,mem=2GB \
      -l walltime=2:00:00:00 \
      -j oe \
      -o /home/nblauch/git/autism/log/${full_log}.log \
      -M blauch@cmu.edu \
      -m ae \
      -q default
    done
  done
done
