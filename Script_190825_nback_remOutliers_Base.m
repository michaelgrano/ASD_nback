load('/Users/Michael/Documents/CMU/Research/ASD/nback/190825_nback_pupilpeakmod_Base.mat')

%delete subjects with no data

autBASE = groups(aut);
conBASE = groups(con);

save('190825_nback_remOutliers_Base.mat','autBASE','conBASE');

%perform same analyses on all groups (autism/control) (distractors/no
%distractors)

function sub = groups(sub)

    base = [sub.BASE];
    sub = outliers(base,sub,6)';

end


%remove outliers

function subw = outliers(DV,sub,f)

    upper = nanmean(DV)+3*nanstd(DV); %upper limit for outlier detection
    lower = nanmean(DV)-3*nanstd(DV); %lower limit for outlier detection
    subCell = struct2cell(sub); %data in cell form
    fieldNames = fieldnames(sub); %fieldnames of structure with data 
    for r = 1:size(subCell,3)
        
        %for each data point
        if subCell{f,1,r} < lower
            %if the value is too low
            subCell{f,1,r} = NaN;
        end
        if subCell{f,1,r} > upper
            %if the value is too high
            subCell{f,1,r} = NaN;
        end
        
    end
    subw = cell2struct(subCell,fieldNames,1); %winsorised data
    subw = subw';

end