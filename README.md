Code to run experiment:
nback_Experiment.m

Code to import/organize data into MATLAB from .csv files:
nback_LoadFromExcel.m (this loads data from the experiment task)
nback_LoadFromExcel_Base.m (this loads data from the baseline recordings prior to starting the experiment)

Code to preprocess data, do behavioral analysis, generate IRFs:
nback_pupilpeakmod.m (this is for the experiment task)
nback_pupilpeakmod_Base.m (this is for the baseline recordings)

Code to remove outliers:
nback_remOutliers.m (this is for the experiment task)
nback_remOutliers_Base.m (this is for the baseline recordings)

Code to summarize preprocessed data and primary dependent measures of interest in a .csv file:
nback_export.CSV.m

Statistical analyses:
nback_Analysis.ipynb
