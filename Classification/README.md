
# Classification analyses for ASD nback experiment

### Requirements:
- `python 3.6+`
  - `sklearn`
  - `matplotlib`
  - `numpy`
  - `scipy`
  - `pickle`
  - `tqdm`
  - `argparse`
  - `os`

### We have pre-generated the data used in classification analyses. To get these `.mat` files directly from `data/nback_remOutliers.mat` , do:
- ` matlab setup_data.m `

### To run classification analyses on the median absolute deviation, use `classify_MAD.slurm`

1)
- If you have access to SLURM-based cluster computing, clone this repository to your cluster and edit the header information in run_slurm.sbatch accordingly.
- If you do not have access to a cluster, change line 9 of `classify_MAD.slurm` to `use_slurm=0`
2) ` . classify_MAD.slurm`

Alternatively, you can call the python script directly with a command line interface (do `python classify_MAD.py --help` to see options)
