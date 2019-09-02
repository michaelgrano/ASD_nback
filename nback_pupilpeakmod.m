%preallocation

subjects = [22 23 24 25 26 27 28 29 30 31 32 33 34 36 37 40 41 43 44 45 46 48 49 50 51 52 53 54 56 57 58 59 60 61 63 64 65 66 68 69 70 72 73 74 75 76 77]; %subjects being ran thru code

for sense = 0:1
    
    %for non-distractor and distractor trials
    folder = '/Users/Michael/Documents/CMU/Research/ASD/nback/'; %folder where data is located
    downsampleRate = 25; %downsample rate
    putativeIRFlength = 4; %presumed IRF length (s)
    samplerate = 1000; %sample rate of EyeLink (Hz)
    filtercutoff = 50; %Butterworth filter cut-off
    butterorder = 2; %Butterworth filter order

    %get going running subjects thru code

    for n = subjects(1:length(subjects))

        %for each participant

        %load and preprocess data

        tic
        if sense == 0
            %if there were no stimuli
            load(strcat(folder,'MatDataSub',num2str(n),'190818')); %load EyeLink data
        else
            %if there were stimuli
            load(strcat(folder,'StimMatDataSub',num2str(n),'190818'));%load EyeLink data
        end
        sub = thedata; %structure with all of the EyeLink data
        
        %data check and figuring out where each block starts/ends
        
        numGaps = 0; %counter for number of gaps in timestamps
        blockStart = [1 NaN NaN]; %start index for each block
        blockEnd = [NaN NaN length(sub.timei)]; %end index for each block
        for r = 1:length(sub.timei)-1
            
            %for each time point
            if sub.timei(r)+1 ~= sub.timei(r+1)
                %if there is a gap in the timestamps
                if numGaps < 2
                    %if this is a normal number of gaps because there is a
                    %total of 3 blocks
                    numGaps = numGaps+1;
                    blockStart(numGaps+1) = r+1;
                    blockEnd(numGaps) = r;
                else
                    disp('There is a boo boo with the timestamps!'); %display error message
                    r
                    return
                end
            end
            
        end
        
        %pupillometry analysis    
        
        sub = myBlink_interpolate(sub); %interpolate out blinks: this uses the Eyelink detected blinks and finds other blinks that the firmware may not have picked up and interpolates over them

        meanpsize1 = mean(sub.psize_all(blockStart(1):blockEnd(1))); %mean pupil size of first block
        meanpsize2 = mean(sub.psize_all(blockStart(2):blockEnd(2))); %mean pupil size of second block
        meanpsize3 = mean(sub.psize_all(blockStart(3):blockEnd(3))); %mean pupil size of third block
        psize1 = (sub.psize_all(blockStart(1):blockEnd(1))-meanpsize1)./meanpsize1;
        psize2 = (sub.psize_all(blockStart(2):blockEnd(2))-meanpsize2)./meanpsize2;
        psize3 = (sub.psize_all(blockStart(3):blockEnd(3))-meanpsize3)./meanpsize3;
        sub.psize_all = [psize1 psize2 psize3]; %normalised pupil sizes (converted to % change relative to mean pupil size per session)
        [b,a] = butter(butterorder,0.03/filtercutoff,'high'); %transfer coefficients for the high-pass filter
        sub.psize_all = [repmat(sub.psize_all(1,1),1,length(sub.psize_all)),sub.psize_all,repmat(sub.psize_all(1,end),1,length(sub.psize_all))]; %convolution boundary handling: pad with first and last sample of signal; important so you don't get onset artifacts, which can create false correlations between frequency components
        sub.psize_all = filtfilt(b,a,sub.psize_all); %perform filtering with filtfilt, to approximately correct for the phase response distortion; this is important because the entire timeseries will shift otherwise, which will mess up your analyses
        sub.psize_all = sub.psize_all(1+length(sub.psize_all)/3:length(sub.psize_all)/3*2); %now you cut off the padding from the boundary handling. 
        sub.psize_all = sub.psize_all';
        psize1 = decimate(sub.psize_all(blockStart(1):blockEnd(1)),downsampleRate); %downsampled pupil data (first block)
        psize2 = decimate(sub.psize_all(blockStart(2):blockEnd(2)),downsampleRate); %downsampled pupil data (second block)
        psize3 = decimate(sub.psize_all(blockStart(3):blockEnd(3)),downsampleRate); %downsampled pupil data (third block)
        sub.psize_all = [psize1; psize2; psize3]; %downsampled pupil sizes
        sub.trialTimeVector_hit = zeros(1,length(sub.psize_all)); %vector where there is a 1 at the stimulus presentation preceding a hit
        sub.trialTimeVector_FA = zeros(1,length(sub.psize_all)); %vector where there is a 1 at the stimulus presentation preceding a FA (estimated @ 1 s)
        sub.trialTimeVector_miss = zeros(1,length(sub.psize_all)); %vector where there is a 1 at the stimulus presentation preceding a miss
        sub.trialTimeVector_aud = zeros(1,length(sub.psize_all)); %vector where there is a 1 at the onset of the distractor stimulus presentation
        
        sub.numRepeats = 0; %number of letter repeats
        sub.numHits = 0; %number of hits
        sub.numFAs = 0; %number of FAs
        sub.numMisses = 0; %number of misses
        outOfTime = 0; %number of letter repeats where there was not enough time to make a response
        blockCount = 1; %counter for what block is being processed
        for r = 1:length(sub.time_all)

            %for each time point
            if r > blockEnd(blockCount)
                %if block has ended
                blockCount = blockCount+1;
            end
            tempstr = sub.time_all{r};
            if strfind(tempstr,'Letter repeated') > 0
                %if there is a letter repeat
                sub.numRepeats = sub.numRepeats+1;
                miss = 1; %switch to determine if a miss has occurred
                if r+putativeIRFlength*1000-1 <= blockEnd(blockCount)
                    %if this isn't right at the end of the block
                    for rr = r:length(sub.time_all)

                        %for each subsequent time point
                        tempstr2 = sub.time_all{rr};
                        if strfind(tempstr2,'Correct') > 0
                            %if a correct response was made
                            sub.trialTimeVector_hit(round(r/downsampleRate)) = 1;
                            sub.numHits = sub.numHits+1;
                            if strfind(tempstr2,'Letter repeated') > 0
                                %if the button press was just about
                                %simultaneous with the letter repeat
                                sub.RTind(sub.numHits,1) = 0;                                
                            else
                                sub.RTind(sub.numHits,1) = (str2double(regexp(tempstr2,'\d+','match'))-str2double(regexp(tempstr,'\d+','match')))/1000; %RTs for individual trials
                            end
                            miss = 0;
                            break
                        elseif rr > r
                            %if the button press isn't simultaneous with
                            %the letter repeat
                            if strfind(tempstr2,'Letter repeated') > 0
                                %if the next letter has repeated without a
                                %button press (i.e., a miss)
                                break
                            end
                        end

                    end
                else
                    outOfTime = outOfTime+1;
                    if blockCount == 3
                        %if this is the last block
                        break
                    end
                    miss = 0;
                end
                if miss == 1
                    %if there was a miss
                    sub.trialTimeVector_miss(round(r/downsampleRate)) = 1;
                    sub.numMisses = sub.numMisses+1;
                end
            elseif strfind(tempstr,'Wrong') > 0
                %if there is a FA
                sub.numFAs = sub.numFAs+1;
                if r > blockStart(blockCount)+999
                    %if the FA is not too close to the start of the block
                    if r+putativeIRFlength*1000-1 <= blockEnd(blockCount)
                        %if this isn't right at the end of the block
                        sub.trialTimeVector_FA(round((r-1000)/downsampleRate)) = 1;
                    end
                else
                    sub.trialTimeVector_FA(blockStart(blockCount)) = 1;
                end
            elseif strfind(tempstr,'No test') > 0
                %if this is the onset of the distractor stimuli
                if r+putativeIRFlength*1000-1 <= blockEnd(blockCount)
                    %if this isn't right at the end of the block
                    if round(r/downsampleRate) > 0
                        %if this is not too close to the start of the block
                        sub.trialTimeVector_aud(round(r/downsampleRate)) = 1;
                    else
                        sub.trialTimeVector_aud(1) = 1;
                    end
                end
            end

        end
        
        %behavioral analysis
        
        sub.numRepeats = sub.numRepeats-outOfTime;
        sub.hitRate = sub.numHits./sub.numRepeats; %hit rate
        if sense == 0
            %no distractors
            sub.faRate = sub.numFAs./(527+529+528-sub.numRepeats); %FA rate
        else
            %distractors
            sub.faRate = sub.numFAs./(533+534+540-sub.numRepeats);
        end
        sub.dprime = norminv(sub.hitRate)-norminv(sub.faRate); %d'
        sub.criterion = -(norminv(sub.hitRate)+norminv(sub.faRate))/2; %C
        sub.RT = mean(sub.RTind); %mean RT
        
        %pupillometry deconvolution etc

        widthDN = putativeIRFlength.*samplerate./downsampleRate; %in samples, downsampled from eyetracker sampling rate
        len_pres = length(sub.trialTimeVector_hit); %length of trialTimeVector_pres
        for ii = 1:widthDN

            %for each column of the design matrix
            cmatrix_hit(1:len_pres,ii) = [zeros(1,ii-1) sub.trialTimeVector_hit(1:len_pres-ii+1)]'; %design matrix (hits)
            cmatrix_FA(1:len_pres,ii) = [zeros(1,ii-1) sub.trialTimeVector_FA(1:len_pres-ii+1)]'; %design matrix (FAs)
            cmatrix_miss(1:len_pres,ii) = [zeros(1,ii-1) sub.trialTimeVector_miss(1:len_pres-ii+1)]'; %design matrix (misses)
            cmatrix_aud(1:len_pres,ii) = [zeros(1,ii-1) sub.trialTimeVector_aud(1:len_pres-ii+1)]'; %design matrix (aud)

        end
        cmatrix_hit = cmatrix_hit(1:length(sub.psize_all),:);
        cmatrix_FA = cmatrix_FA(1:length(sub.psize_all),:);
        cmatrix_miss = cmatrix_miss(1:length(sub.psize_all),:);
        cmatrix_aud = cmatrix_aud(1:length(sub.psize_all),:);
        cmatrix_all = horzcat(cmatrix_hit,cmatrix_FA,cmatrix_miss); %combined design matrix
        sub.deconvolved = sub.psize_all'*pinv(cmatrix_all)'; %deconvolved response (first 1/3 is for hits, second 1/3 is for FAs, third 1/3 is for misses)
        sub.taskIrf_hit = sub.deconvolved(1,1:1000/downsampleRate*putativeIRFlength); %deconvolved response (hits)
        sub.taskIrf_FA = sub.deconvolved(1,1000/downsampleRate*putativeIRFlength+1:2*1000/downsampleRate*putativeIRFlength); %deconvolved response (FAs)
        sub.taskIrf_miss = sub.deconvolved(1,2*1000/downsampleRate*putativeIRFlength+1:3*1000/downsampleRate*putativeIRFlength); %deconvolved response (misses)
        sub.taskIrf_aud = sub.psize_all'*pinv(cmatrix_aud)'; %deconvolved response (distractor stimuli)

        [sub.taskIrf_hit,sub.taskIrf_hit_mod,sub.taskIrf_hit_max,sub.taskIrf_hit_abs] = modIrf(sub.taskIrf_hit,sub.trialTimeVector_hit,putativeIRFlength,downsampleRate); %analyse the IRF (hits)
        [sub.taskIrf_FA,sub.taskIrf_FA_mod,sub.taskIrf_FA_max,sub.taskIrf_FA_abs] = modIrf(sub.taskIrf_FA,sub.trialTimeVector_FA,putativeIRFlength,downsampleRate); %analyse the IRF (FAs)
        [sub.taskIrf_miss,sub.taskIrf_miss_mod,sub.taskIrf_miss_max,sub.taskIrf_miss_abs] = modIrf(sub.taskIrf_miss,sub.trialTimeVector_miss,putativeIRFlength,downsampleRate); %analyse the IRF (misses)
        [sub.taskIrf_aud,sub.taskIrf_aud_mod,sub.taskIrf_aud_max,sub.taskIrf_aud_abs] = modIrf(sub.taskIrf_aud,sub.trialTimeVector_aud,putativeIRFlength,downsampleRate); %analyse the IRF (distractor stimuli)
        
        if sense == 0
            %if there are no stimuli
            sub_analysed(n) = sub; %output from this code
        else
            sub_analysed_a(n) = sub; %output from this code (stimulus trials)
        end
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
        if x == 1 && sense == 0
            %if subject does have autism and if no stim trials
            aut(a) = sub_analysed(s); %code output (autism only)
            a = a+1;
        elseif x == 1 && sense == 1
            %if subject does have autism and if stim trials
            aut_a(a) = sub_analysed_a(s); %code output (autism only) (distractor trials)
            a = a+1;
        elseif x == 0 && sense == 0
            %if subject is a control and if no stim strials
            con(c) = sub_analysed(s); %code output (control only)
            c = c+1;
        else
            %if subject is a control and if stim trials
            con_a(c) = sub_analysed_a(s); %code output (control only) (distractor trials)
            c = c+1;
        end

    end
    
end

aut = rmfield(aut,{'time_all','psize_all'});
aut_a = rmfield(aut_a,{'time_all','psize_all'});
con_a = rmfield(con_a,{'time_all','psize_all'});
con = rmfield(con,{'time_all','psize_all'});

save('190824_nback_pupilpeakmod.mat','aut','con','aut_a','con_a');

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

%perform same thing for all types of IRFs

function [taskIrf,taskIrf_mod,taskIrf_max,taskIrf_abs] = modIrf(taskIrf,trialTimeVector,putativeIRFlength,downsampleRate) %analyse the IRF (presentations)

    taskIrf_mod = mad(taskIrf,1); %mad of deconvolved response
    taskIrf_max = max(taskIrf); %max of deconvolved response
    taskIrf_abs = max(taskIrf)-min(taskIrf); %absolute difference of max - min of deconvolved response
    if sum(trialTimeVector) == 0
        %if there were not enough events to generate an IRF
        taskIrf = nan(1,putativeIRFlength*1000/downsampleRate);
        taskIrf_mod = NaN;
        taskIrf_max = NaN;
        taskIrf_abs = NaN;
    end

end
