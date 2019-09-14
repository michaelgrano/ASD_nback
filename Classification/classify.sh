
source activate mne
cd ~/git/autism

nd='--compute-nulldist'
# irf='--short-irf'
irf=''
# nd=''
# binspace=$1
# outer_cv=$2

python classify.py --bin-width $binwidth --bin-space $binspace --outer-cv $outer_cv --dtype $dtype --overwrite $nd $irf $alt_tag
