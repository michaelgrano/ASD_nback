load('/Users/Michael/Documents/CMU/Research/ASD/nback/190824_nback_pupilpeakmod.mat')
putativeIRFlength = 4;
downsampleRate = 25;

aut = groups(aut,putativeIRFlength,downsampleRate);
aut_a = groups(aut_a,putativeIRFlength,downsampleRate);
con = groups(con,putativeIRFlength,downsampleRate);
con_a = groups(con_a,putativeIRFlength,downsampleRate);

save('190825_nback_remOutliers.mat','aut','con','aut_a','con_a');

%perform same analyses on all groups (autism/control) (distractors/no
%distractors)

function sub = groups(sub,putativeIRFlength,downsampleRate)

    hit = [sub.taskIrf_hit_mod];
    FA = [sub.taskIrf_FA_mod];
    miss = [sub.taskIrf_miss_mod];
    aud = [sub.taskIrf_aud_mod];
    sub = outliers(hit,sub,25,putativeIRFlength,downsampleRate)';
    sub = outliers(FA,sub,28,putativeIRFlength,downsampleRate)';
    sub = outliers(miss,sub,31,putativeIRFlength,downsampleRate)';
    sub = outliers(aud,sub,34,putativeIRFlength,downsampleRate)';
    Dprime = [sub.dprime];
    Criterion = [sub.criterion];
    maxD = max(Dprime(~isinf(Dprime)));
    maxC = max(Criterion(~isinf(Criterion)));
    minC = min(Criterion(~isinf(Criterion)));
    sub = blanks(sub,17,maxD)';
    sub = blanks(sub,18,maxC,minC)';
    sub = rmfield(sub, {'blink_start','blink_end','timei','trialTimeVector_hit','trialTimeVector_FA','trialTimeVector_miss','trialTimeVector_aud','numHits','numFAs','numMisses','numRepeats','deconvolved'});

end


%remove outliers

function subw = outliers(DV,sub,f,putativeIRFlength,downsampleRate)

    upper = nanmean(DV)+3*nanstd(DV); %upper limit for outlier detection
    lower = nanmean(DV)-3*nanstd(DV); %lower limit for outlier detection
    subCell = struct2cell(sub); %data in cell form
    fieldNames = fieldnames(sub); %fieldnames of structure with data 
    for r = 1:size(subCell,3)
        
        %for each data point
        if subCell{f,1,r} < lower
            %if the value is too low
            subCell{f/3+38/3,1,r} = nan(1,putativeIRFlength*1000/downsampleRate);
            subCell{f,1,r} = NaN;
        end
        if subCell{f,1,r} > upper
            %if the value is too high
            subCell{f/3+38/3,1,r} = nan(1,putativeIRFlength*1000/downsampleRate);
            subCell{f,1,r} = NaN;
        end
        
    end
    subw = cell2struct(subCell,fieldNames,1); %winsorised data
    subw = subw';

end


%remove blanks

function subw = blanks(sub,f,max,min)

    subCell = struct2cell(sub);
    fieldNames = fieldnames(sub); %fieldnames of structure with data 
    for r = 1:size(subCell,3)
        
        %for each data point
        if subCell{f,1,r} == Inf
            %if the data point is infinity
            subCell{f,1,r} = max;
        elseif subCell{f,1,r} == -Inf
            %if the data point is -infinity
            subCell{f,1,r} = min;
        end
        
    end
    subw = cell2struct(subCell,fieldNames,1);
    subw = subw';

end
