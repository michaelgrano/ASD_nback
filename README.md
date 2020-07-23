Code to run experiment:

      nback_Experiment.m


Code to import/organize data into MATLAB from .csv files:

      nback_LoadFromExcel.m (this loads data from the experiment task)

      nback_LoadFromExcel_Base.m (this loads data from the baseline recordings prior to starting the experiment)


Code to preprocess data, do behavioral analysis, generate IRFs:

      nback_pupilpeakmod.m (this is for the experiment task) (_IRF4 version uses transformed RTs)

      nback_pupilpeakmod_Base.m (this is for the baseline recordings)
      
      Blink interpolation code has been adapted from the following source:

% This code reproduces the analyses in the paper
% Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven 
% by decision uncertainty and alters serial choice bias. 
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai, 2016
% anne.urai@gmail.com


Code to remove outliers:

      nback_remOutliers.m (this is for the experiment task)

      nback_remOutliers_Base.m (this is for the baseline recordings)


Code to summarize preprocessed data and primary dependent measures of interest in a .csv file:

      nback_export.CSV.m


Code to compute ratios of pupil response amplitudes between conditions:

      nback_Ratios.m
      
      
Output of nback_export.CSV.m can be imported into nback_Analysis_IRF4_Final.ipynb to replicate primary results of analyses in paper.

Instructions for classification code is in Classification folder. Use ClassificationUSETHISCODE for code and data.
