
source activate mne
cd ~/git/autism

nd='--compute-nulldist'
# irf='--short-irf'
irf=''

python classify_MAD.py --dtype $dtype $alt_tag $nd --outer-cv $outer_cv --measure $mad_type $overwrite $irf
