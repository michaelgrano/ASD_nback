spreadsheet = '/Users/Michael/Documents/CMU/Research/ASD/nback/190825_SubjectData'; %name of spreadsheet where subject info. is stored

[~,~,subjectNumber] = xlsread(spreadsheet,1,'C2:C48'); %subject number
[~,~,diagnosis] = xlsread(spreadsheet,1,'E2:E48'); %diagnosis
[~,~,age] = xlsread(spreadsheet,1,'F2:F48'); %age
[~,~,female] = xlsread(spreadsheet,1,'G2:G48'); %gender
[~,~,handedness] = xlsread(spreadsheet,1,'I2:I48'); %handedness?
[~,~,caffeine] = xlsread(spreadsheet,1,'K2:K48'); %caffeine?
[~,~,meds] = xlsread(spreadsheet,1,'L2:L48'); %meds that affect the adrenergic system?
[~,~,glasses] = xlsread(spreadsheet,1,'M2:M48'); %glasses?
[~,~,pupilContrast] = xlsread(spreadsheet,1,'N2:N48'); %pupil contrast?
[~,~,corneaContrast] = xlsread(spreadsheet,1,'O2:O48'); %cornea contrast?
[~,~,ADOSComm] = xlsread(spreadsheet,1,'U2:U48'); %ADOSComm
[~,~,ADOSSoc] = xlsread(spreadsheet,1,'V2:V48'); %ADOSSoc
[~,~,ADOSSocComm] = xlsread(spreadsheet,1,'W2:W48'); %ADOSSocComm
[~,~,ADOSBeh] = xlsread(spreadsheet,1,'X2:X48'); %ADOSBeh
[~,~,VIQ] = xlsread(spreadsheet,1,'Y2:Y48'); %VIQ
[~,~,PIQ] = xlsread(spreadsheet,1,'Z2:Z48'); %PIQ
[~,~,FIQ] = xlsread(spreadsheet,1,'AA2:AA48'); %FIQ

subjectNumber = cell2mat(subjectNumber); %convert cell into array
diagnosis = cell2mat(diagnosis); %convert cell into array
age = cell2mat(age); %convert cell into array
female = cell2mat(female); %convert cell into array
handedness = cell2mat(handedness); %convert cell into array
caffeine = cell2mat(caffeine); %convert cell into array
meds = cell2mat(meds); %convert cell into array
glasses = cell2mat(glasses); %convert cell into array
pupilContrast = cell2mat(pupilContrast); %convert cell into array
corneaContrast = cell2mat(corneaContrast); %convert cell into array
ADOSComm = cell2mat(ADOSComm); %convert cell into array
ADOSSoc = cell2mat(ADOSSoc); %convert cell into array
ADOSSocComm = cell2mat(ADOSSocComm); %convert cell into array
ADOSBeh = cell2mat(ADOSBeh); %convert cell into array
VIQ = cell2mat(VIQ); %convert cell into array
PIQ = cell2mat(PIQ); %convert cell into array
FIQ = cell2mat(FIQ); %convert cell into array

load('/Users/Michael/Documents/CMU/Research/ASD/nback/190825_nback_remOutliers.mat')
load('/Users/Michael/Documents/CMU/Research/ASD/nback/190825_nback_remOutliers_Base.mat')

subjectNumber = repeatInfo(subjectNumber);
diagnosis = repeatInfo(diagnosis);
age = repeatInfo(age);
female = repeatInfo(female);
handedness = repeatInfo(handedness);
caffeine = repeatInfo(caffeine);
meds = repeatInfo(meds);
glasses = repeatInfo(glasses);
pupilContrast = repeatInfo(pupilContrast);
corneaContrast = repeatInfo(corneaContrast);
ADOSComm = repeatInfo(ADOSComm);
ADOSSoc = repeatInfo(ADOSSoc);
ADOSSocComm = repeatInfo(ADOSSocComm);
ADOSBeh = repeatInfo(ADOSBeh);
VIQ = repeatInfo(VIQ);
PIQ = repeatInfo(PIQ);
FIQ = repeatInfo(FIQ);
baseline = [[autBASE.BASE]'; [conBASE.BASE]'; [autBASE.BASE]'; [conBASE.BASE]'; [autBASE.BASE]'; [conBASE.BASE]'; [autBASE.BASE]'; [conBASE.BASE]'; [autBASE.BASE]'; [conBASE.BASE]'; [autBASE.BASE]'; [conBASE.BASE]'; [autBASE.BASE]'; [conBASE.BASE]'];
dprime = [[aut.dprime]'; [con.dprime]'; [aut.dprime]'; [con.dprime]'; [aut.dprime]'; [con.dprime]'; [aut_a.dprime]'; [con_a.dprime]'; [aut_a.dprime]'; [con_a.dprime]'; [aut_a.dprime]'; [con_a.dprime]'; [aut_a.dprime]'; [con_a.dprime]'];
criterion = [[aut.criterion]'; [con.criterion]'; [aut.criterion]'; [con.criterion]'; [aut.criterion]'; [con.criterion]'; [aut_a.criterion]'; [con_a.criterion]'; [aut_a.criterion]'; [con_a.criterion]'; [aut_a.criterion]'; [con_a.criterion]'; [aut_a.criterion]'; [con_a.criterion]'];
RT = [[aut.RT]'; [con.RT]'; [aut.RT]'; [con.RT]'; [aut.RT]'; [con.RT]'; [aut_a.RT]'; [con_a.RT]'; [aut_a.RT]'; [con_a.RT]'; [aut_a.RT]'; [con_a.RT]'; [aut_a.RT]'; [con_a.RT]'];
peak = [[aut.taskIrf_hit_mod]'; [con.taskIrf_hit_mod]'; [aut.taskIrf_FA_mod]'; [con.taskIrf_FA_mod]'; [aut.taskIrf_miss_mod]'; [con.taskIrf_miss_mod]'; [aut_a.taskIrf_hit_mod]'; [con_a.taskIrf_hit_mod]'; [aut_a.taskIrf_FA_mod]'; [con_a.taskIrf_FA_mod]'; [aut_a.taskIrf_miss_mod]'; [con_a.taskIrf_miss_mod]'; [aut_a.taskIrf_aud_mod]'; [con_a.taskIrf_aud_mod]']; %peak pupil size
hit = [ones(47,1); zeros(94,1); ones(47,1); zeros(141,1)]; %hit?
FA = [zeros(47,1); ones(47,1); zeros(94,1); ones(47,1); zeros(94,1)]; %FA?
miss = [zeros(94,1); ones(47,1); zeros(94,1); ones(47,1); zeros(47,1)]; %miss?
aud = [zeros(282,1); ones(47,1)]; %tones?
distractors = [zeros(47*3,1); ones(47*4,1)]; %condition

tidyData = table(subjectNumber,diagnosis,age,female,handedness,caffeine,meds,glasses,pupilContrast,corneaContrast,ADOSComm,ADOSSoc,ADOSSocComm,ADOSBeh,VIQ,PIQ,FIQ,baseline,dprime,criterion,RT,peak,hit,FA,miss,aud,distractors);
writetable(tidyData,'190825_tidyData.csv')

function var = repeatInfo(var)

    var = [var; var; var; var; var; var; var];

end