%preallocation

subjects = [22 23 24 25 26 27 28 29 30 31 32 33 34 36 37 40 41 43 44 45 46 48 49 50 51 52 53 54 56 57 58 59 60 61 63 64 65 66 68 69 70 72 73 74 75 76 77]; %subjects being ran thru code

folder = '/Users/Michael/Documents/CMU/Research/ASD/nback/'; %folder where data is located
samplerate = 1000; %sample rate of EyeLink (Hz)

%get going running subjects thru code

for n = subjects(1:length(subjects))

    %for each participant

    %load and preprocess data

    tic
    load(strcat(folder,'MatBaseDataSub',num2str(n),'190815')); %load EyeLink data
    sub = thedata; %structure with all of the EyeLink data
    load(strcat(folder,'StimMatBaseDataSub',num2str(n),'190815'));%load EyeLink data
    sub_a = thedata; %structure with all of the EyeLink data (distractors condition)
    sub.timei = [sub.timei; sub_a.timei];
    sub.psize_all = [sub.psize_all; sub_a.psize_all];

    %data check and figuring out where each block starts/ends

    numGaps = 0; %counter for number of gaps in timestamps
    for r = 1:length(sub.timei)-1

        %for each time point
        if sub.timei(r)+1 ~= sub.timei(r+1)
            %if there is a gap in the timestamps
            if numGaps < 5
                %if this is a normal number of gaps because there is a
                %total of 6 blocks
                numGaps = numGaps+1;
            else
                disp('There is a boo boo with the timestamps!'); %display error message
                r
                return
            end
        end

    end

    %pupillometry analysis    

    sub = myBlink_interpolate(sub); %interpolate out blinks: this uses the Eyelink detected blinks and finds other blinks that the firmware may not have picked up and interpolates over them

    sub.BASE = nanmedian(sub.psize_all); %median baseline size
    sub_analysed(n) = sub;
    
    n %subject number
    toc %time to run participant's data thru code

end

autistics = [22 23 24 25 26 28 34 40 41 45 49 50 51 52 54 60 63 64 65 66 68 69 70]; %subject IDs of autistic participants
controls = [27 29 30 31 32 33 35 36 37 43 44 46 48 53 56 57 58 59 61 72 73 74 75 76 77]; %subject IDs of controls
a = 1; %counter for number of subjects with autism
c = 1; %counter for number of controls
for s = subjects(1:length(subjects))

    %for each subject
    x = ismember(s,autistics); %switch to see if subject has autism 
    if x == 1
        %if subject does have autism and if no stim trials
        aut(a) = sub_analysed(s); %code output (autism only)
        a = a+1;
    elseif x == 0
        %if subject is a control and if no stim strials
        con(c) = sub_analysed(s); %code output (control only)
        c = c+1;
    end

end
    
aut = rmfield(aut,{'time_all','psize_all'});
con = rmfield(con,{'time_all','psize_all'});

save('190825_nback_pupilpeakmod_Base.mat','aut','con');

%interpolate blinks and missing data
%Anne Urai, 2016, modified by eli and charlie
%modified 2018 by michael

function sub = myBlink_interpolate(sub)

    %get the data we need

    dat.time = sub.timei'; %timestamps
    dat.pupil = sub.psize_all'; %psizes
    padding = 0.150; % how long before and after do we want to pad? (s)
    data.fsample = 1000; %sample rate
    blinksmp = []; %two-column matrix where the values in the first column are the indices of the start of each blink and those in the second are the indices of the end of each blink
    for iBlink = 1:length(sub.blink_start)

        %for each blink
        if isempty(find(dat.time==sub.blink_start(iBlink),1))
            %if this is empty (i.e., if the block started with a blink)
            dat.time = [dat.time(1:find(dat.time==sub.blink_start(iBlink)-1)) sub.blink_start(iBlink) dat.time(find(dat.time==sub.blink_start(iBlink)-1)+1:end)];
            dat.pupil = [dat.pupil(1:find(dat.time==sub.blink_start(iBlink)-1)) NaN dat.pupil(find(dat.time==sub.blink_start(iBlink)-1)+1:end)];
        elseif isempty(find(dat.time==sub.blink_end(iBlink),1))
            %if this is empty (i.e., if the block ended with a blink)
            dat.time = [dat.time(1:find(dat.time==sub.blink_end(iBlink)-1)) sub.blink_end(iBlink) dat.time(find(dat.time==sub.blink_end(iBlink)-1)+1:end)];
            dat.pupil = [dat.pupil(1:find(dat.time==sub.blink_end(iBlink)-1)) NaN dat.pupil(find(dat.time==sub.blink_end(iBlink)-1)+1:end)];
        end
        blinksmp = cat(1,blinksmp,[find(dat.time==sub.blink_start(iBlink)) find(dat.time==sub.blink_end(iBlink))]);

    end

    % initialize output

    if ~isempty(blinksmp)
        %if there are blinks

        %interpolate EL-defined blinks

        %merge 2 blinks into 1 if they are < 250 ms together (coalesce)

        coalesce = 0.250; %minimum time interval between blinks for each blink to be considered distinct (s)
        for b = 1:size(blinksmp,1)-1

            %for each blink (except for the last blink)
            if blinksmp(b+1,1)-blinksmp(b,2) < coalesce*data.fsample
                %if two blinks are too close to one another in time
                blinksmp(b,2) = blinksmp(b+1,2);
                blinksmp(b+1,:) = NaN;
            end

        end

        %remove those duplicates

        blinksmp(isnan(nanmean(blinksmp,2)),:) = [];

        %pad the blinks

        padblinksmp(:,1) = round(blinksmp(:,1)-padding*data.fsample); %blink start indices (padded)
        padblinksmp(:,2) = round(blinksmp(:,2)+padding*data.fsample); %blink end inedices (padded)

        %avoid idx outside range

        if any(padblinksmp(:) < 1)
            %if the index is outside of range
            padblinksmp(padblinksmp < 1) = 1;
        end
        if any(padblinksmp(:) > length(dat.pupil))
            %if the index is outside of range
            padblinksmp(padblinksmp > length(dat.pupil)) = length(dat.pupil);
        end

        %make the pupil NaN at those points

        for b = 1:size(padblinksmp,1)

            %for each blink
            dat.pupil(padblinksmp(b,1):padblinksmp(b,2)) = NaN;

        end

        %also remove outliers

        dat.pupil(dat.pupil < nanmedian(dat.pupil)-3*nanstd(dat.pupil)) = NaN;
        dat.pupil(dat.pupil > nanmedian(dat.pupil)+3*nanstd(dat.pupil)) = NaN;

        %interpolate linearly

        if sum(isnan(dat.pupil)) > (2/3)*length(dat.pupil)
            %if there are more than 2/3 Nans, return control to invoking function
            disp('WHOA. WAYYYY too many blinks!'); %display error message in command window
            return
        else
            dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)),dat.pupil(~isnan(dat.pupil)),find(isnan(dat.pupil)),'linear');
        end
    end

    %also extrapolate ends

    dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)),dat.pupil(~isnan(dat.pupil)),find(isnan(dat.pupil)),'nearest','extrap');

    %interpolate peak-detected blinks

    assert(~any(isnan(dat.pupil)));
    win = hanning(11); %Hann window
    pupildatsmooth = filter2(win.',dat.pupil,'same'); %pupil data with finite impulse response filter applied
    dat.pupildiff = diff(pupildatsmooth)-mean(diff(pupildatsmooth))/std(diff(pupildatsmooth)); %input signal vector to detect peaks from
    [peaks,loc] = findpeaks(abs(dat.pupildiff),'minpeakheight',3*std(dat.pupildiff),'minpeakdistance',0.5*data.fsample); %peak values and their corresponding time points

    if ~isempty(peaks)
        %if there are peaks

        %convert peaks into blinksmp
        
        newblinksmp = nan(length(peaks),2); %peak-detected blinks
        for p = 1:length(peaks)
            
            %for each peak-detected blink
            newblinksmp(p,1) = loc(p)-2*padding*data.fsample; %peak detected will be eye-opening again
            newblinksmp(p,2) = loc(p)+padding*data.fsample;
            
        end

        %merge 2 blinks into 1 if they are < 250 ms together (coalesce)
        coalesce = 0.250;
        for b = 1:size(newblinksmp,1)-1
            
            %for each peak-detected blink
            if newblinksmp(b+1,1)-newblinksmp(b,2) < coalesce*data.fsample
                %if two peak-detected blinks are too close to one another
                %in time
                newblinksmp(b,2) = newblinksmp(b+1,2);
                newblinksmp(b+1,:) = NaN;
            end
            
        end
        
        %remove those duplicates
        
        newblinksmp(isnan(nanmean(newblinksmp,2)),:) = [];

        %make sure none are outside of the data range
        
        newblinksmp(newblinksmp < 1) = 1;
        newblinksmp(newblinksmp > length(dat.pupil)) = length(dat.pupil);

        %make the pupil NaN at those points
        
        for b = 1:size(newblinksmp,1)
            
            %for each peak-detected blink
            dat.pupil(newblinksmp(b,1):newblinksmp(b,2)) = NaN;
            
        end

        %interpolate linearly
        
        dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)),dat.pupil(~isnan(dat.pupil)),find(isnan(dat.pupil)),'linear');
    end

    %remove remaining nans (probably at the end)
    
    dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)),dat.pupil(~isnan(dat.pupil)),find(isnan(dat.pupil)),'nearest','extrap');

    %store the results
    
    sub.timei = dat.time';
    sub.psize_all = dat.pupil;

end