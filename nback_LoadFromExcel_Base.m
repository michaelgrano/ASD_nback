%preallocation

filelocation = '/Users/Michael/Documents/CMU/Research/ASD/nback/Sensory_Data_C_V/Subject'; %location of file
filelocation_a = '/Users/Michael/Documents/CMU/Research/ASD/nback/Sensory_Data_C_A/Subject'; %location of file (stim. trials)

subjects = [22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 40 41 43 44 45 46 48 49 50 51 52 53 54 56 57 58 59 60 61 63 64 65 66 68 69 70 72 73 74 75 76 77]; %subjects being ran thru code

%get going for each subject

for sense = 0:1

    %for each condition (presence vs. absence of distractors)
    for n = subjects(1:length(subjects))

        %for each participant
        tic

        %load and sort data

        if sense == 0
            %non-stimulus trials
            cspreadsheet = strcat(filelocation,num2str(n),'.xlsx'); %name of spreadsheet where data is stored
        else
            %stimulus trials
            cspreadsheet = strcat(filelocation_a,num2str(n),'.xlsx'); %name of spreadsheet where data is stored
        end
        
        [time1,x1,psize1] = moveOver(cspreadsheet,1);
        [time2,x2,psize2] = moveOver(cspreadsheet,2);
        [time3,x3,psize3] = moveOver(cspreadsheet,3);
        sub(n).time_all = [time1; time2; time3]; %all timestamp data
        sub(n).x_all = [x1; x2; x3]; %all x-coordinate data
        sub(n).psize_all = [psize1; psize2; psize3]; %all pupil size data

        %format data

        sub(n).num_blinks = 0; %number of times a participant blinks
        sub(n).num_saccades = 0; %number of times a participant saccades
        rtodelete = []; %rows to delete because the EyeLink adds in an extra row for the same time point in the output file every time that a message appears
        sub(n).blink_start = []; %indices at which blinks start
        sub(n).blink_end = []; %indices at which blinks end
        for r = 1:length(sub(n).time_all)

            %for each time point
            tempstr = sub(n).time_all{r};
            if strfind(tempstr,'SBLINK') > 0
                %if a blink starts
                rtodelete = [rtodelete r];
                if isnumeric(sub(n).x_all{r})
                    %if this is not a text stamp (i.e., not a blink
                    %occurring at the start of the block)
                    sub(n).blink_start = [sub(n).blink_start str2double(regexp(tempstr,'\d+','match'))];
                else
                    sub(n).blink_start = [sub(n).blink_start str2double(regexp(sub(n).x_all{r},'\d+','match'))];
                end
                sub(n).num_blinks = sub(n).num_blinks+1;
            elseif strfind(tempstr,'EBLINK') > 0
                %if a blink ends
                rtodelete = [rtodelete r];
                if isnumeric(sub(n).x_all{r})
                    %if this is not a text stamp (i.e., not a blink
                    %occurring at the end of the block)
                    sub(n).blink_end = [sub(n).blink_end sub(n).x_all{r}];
                else
                    sub(n).blink_end = [sub(n).blink_end str2double(regexp(sub(n).x_all{r},'\d+','match'))];
                end
            elseif  strfind(tempstr,'SSACC') > 0
                %if this is a time point at which a saccade begins
                rtodelete = [rtodelete r];
                sub(n).num_saccades = sub(n).num_saccades+1;
            elseif ischar(tempstr)
                %if this is still text
                if ~strcmp(tempstr,'MSG')
                    %if this is not a message line
                    rtodelete = [rtodelete r];
                end
            end
            if sub(n).psize_all{r} == 0
                %if a blink/artifact occurrs
                sub(n).x_all{r} = NaN;
                sub(n).psize_all{r} = NaN;
            end
            if isnumeric(sub(n).psize_all{r}) == 1
                %if at this time point, a useable pupil size is recorded
                sub(n).psize_all{r} = 2.*(sub(n).psize_all{r}./pi).^0.5; %converts area to diameter
            end

        end
        sub(n).time_all(rtodelete) = [];
        sub(n).psize_all(rtodelete) = [];
        sub(n).x_all(rtodelete) = [];
        timei = sub(n).time_all; %timestamps WITHOUT messages
        rtodelete = [];
        for r = 1:length(sub(n).x_all)

            %for each time point
            tempstr = sub(n).x_all{r};
            isMSG = sub(n).time_all{r}; %should be string 'MSG', unless one message immediately preceded another and that message was put here, and this line should not be deleted
            if ischar(tempstr)
                %if a message is here
                if strfind(isMSG,'MSG') > 0
                    %if this line doesn't have relevant info already
                    %displaced in the first column
                    sub(n).time_all{r+1} = sub(n).x_all{r};
                    if r+1 > length(timei)
                        %if this is the last line
                        sub(n).x_all{r+1} = NaN;
                        sub(n).psize_all{r+1} = NaN;
                        timei{r+1} = str2double(regexp(sub(n).x_all{r},'\d+','match'));
                    end
                else
                    sub(n).time_all{r+1} = strcat(sub(n).time_all{r},sub(n).x_all{r});
                end
                rtodelete = [rtodelete r];
            end

        end
        sub = rmfield(sub,{'x_all'});
        sub(n).time_all(rtodelete) = [];
        sub(n).psize_all(rtodelete) = [];
        timei(rtodelete) = [];
        sub(n).psize_all = cell2mat(sub(n).psize_all); %converted to matrix
        sub(n).timei = cell2mat(timei); %converted to matrix

        %save and finish

        thedata = sub(n); %data to be saved
        if sense == 0
            %if there are no stimuli
            save(strcat('MatBaseDataSub',num2str(n),'190815.mat'),'thedata'); %save data
        else
            %if there are stimuli
            save(strcat('StimMatBaseDataSub',num2str(n),'190815.mat'),'thedata');
        end

        n %subject number
        toc %time to run participant's data thru code

    end
    
end

%move data over from Excel into MATLAB, from each Excel sheet

function [time,x,psize] = moveOver(cspreadsheet,sheet)

    [~,~,time_all] = xlsread(cspreadsheet,sheet,'A:A'); %column 1 of data (timestamps and messages)
    [~,~,x_all] = xlsread(cspreadsheet,sheet,'B:B'); %column 2 of data (x-coords and messages)
    [~,~,psize_all] = xlsread(cspreadsheet,sheet,'D:D'); %column 4 of data (psize)
    blinkExists1 = 0; %switch for existence of a blink start
    blinkExists2 = 0; %switch for existence of a blink end
    for r = 1:length(x_all)
        
        %for each time point
        tempstr = x_all{r}; %string if EyeLink sent a trigger message at this time point
        if strfind(tempstr,'Baseline Recording Start') > 0
            %if the block begins here
            design_start = r; %row index at which the block starts
            for rr = r+1:length(x_all)
                
                %for each time point subsequent to the start of the block
                tempstr = x_all{rr};
                if strfind(tempstr,'Baseline Recording End') > 0
                    design_end = rr; %row index at which the block ends
                    break
                end
                
            end
            endloop = 0; %switch to end the loop that is about to commence
            for rr = r+1:design_end
                
                %for each time point subsequent to the start of the block
                tempstr2 = time_all{rr}; %string if EyeLink sent a trigger message at this time point
                if strfind(tempstr2,'SBLINK') > 0
                    %if a blink is started
                    blinkExists1 = 1;
                    tempsblink = rr; %index at which the first blink is started once the block has started
                    for rrr = r+1:design_end
                        
                        %for each time point subsequent to the start of the
                        %block
                        tempstr3 = time_all{rrr}; %string if EyeLink sent a trigger message at this time point
                        if strfind(tempstr3,'EBLINK') > 0
                            %if a blink is ended
                            blinkExists2 = 1;
                            tempeblink = rrr; %index at which the first blink is ended following the start of the block
                            if tempeblink < tempsblink
                                %if the subject started a blink before the
                                %start of the block but did not complete
                                %the blink until after the start of the
                                %block
                                time_all{r} = 'SBLINK';
                            end
                            endloop = 1;
                            break
                        end
                        
                    end
                end
                if endloop == 1
                    %if this loop needs to end
                    break
                end
                
            end
            break
        end
        
    end
    if blinkExists1 == 1
        %if a blink was started during the block
        for r = design_end:-1:design_start

            %for each time point starting at the end of the block
            tempstr = time_all{r};
            if strfind(tempstr,'SBLINK') > 0
                %if a blink is started
                tempsblink2 = r; %index at which the last blink is started in a block
                break
            end

        end
    else
        tempsblink2 = 1;
    end
    if blinkExists2 == 1
        %if a blink was ended during the block
        for r = design_end:-1:design_start

            %for each time point starting at the end of the block
            tempstr = time_all{r};
            if strfind(tempstr,'EBLINK') > 0
                %if a blink is ended
                tempeblink2 = r; %index at which the last blink is ended in a block
                break
            end

        end
    elseif blinkExists1 == 1
        %if a blink was started during the block
        tempeblink2 = 0;
    else
        tempeblink2 = 2;
    end
    if tempsblink2 > tempeblink2
        %if a blink is started at the end of a block but not finished
        time_all{design_end} = 'EBLINK';
    end
    time = time_all(design_start:design_end,1); %truncated column 1 data to start of block
    x = x_all(design_start:design_end,1); %truncated column 2 data to start of block
    psize = psize_all(design_start:design_end,1); %truncated column 4 data to start of block

end
