% 08/02/10
% Auditory adaptation experiment.
% Play two consecutive tones to both ears simultaneously.
% Each trial contains an adapter with five beeps and then a tester with one
% beep followed by the iti. Subjects perform an orthogonal visual letter
% one back task at fixation.

% Check that num_trials = 36    (12 seconds per trial on average)
% Use in magnet with TR = 1.5 sec
% Total run time is: 15(blank) + 432(trials) + 9(blank) = 456sec = 304TRs

%07/2019
%Adapted to use code to assess pupillometry on blocks in the presence and
%absence of distractor stimuli

Screen('Preference','SkipSyncTests',1); %Skip Psychtoolbox sync testing


clear all %clear all
KbName('UnifyKeyNames') %important to get keys for Windows
sense = 1; %switch for block type: set this to 0 for no distractors, 1 for distractors
eyetracking = 1; %switch for eyetracking
path = [cd '/']; %path
sub = 'test'; %VARIABLE NOT USED FOR CURRENT PURPOSES

number = input('Subject number? '); %subject ID

AssertOpenGL; %Initialize
num_trials = 24;     % Number of trials in the experiment (code set for 35). %2016: for current purposes, this is the number of auditory tone onsets in the distractors condition
freq = 44100;        % Frequency at which to present the stimulus
freq_t1 = 600;      % Frequency of adapter tone
freq_t2 = 400;       % Frequency of tester tone %VARIABLE NOT USED FOR CURRENT PURPOSES
dt = 1/freq;         % Step between each sample - for tones.
letter_length = 0.5;  % Time to show letter at fixation in seconds
beep_length = 0.15;    % Beep length in seconds
response_time = 1;  % Allowed response time in seconds
mon_width   = 38;   % horizontal dimension of viewable screen (cm)
v_dist      = 60;   % viewing distance (cm)
fix_r       = 0.7;  % radius of fixation point (deg)

%make distractor tones

t_beep = [dt:dt:beep_length]; % Beep length.
t1 = sin(2*pi*freq_t1*t_beep);
t2 = sin(2*pi*freq_t2*t_beep);

% cosine ramp
Tattack = 0.005;
A=(0:dt:Tattack)/Tattack;
Tfade=(pi/(length(A)-.5));
RaisedCosine=cos(pi:Tfade:3*pi)+1;
RaisedCosineNormSquare=(RaisedCosine/max(RaisedCosine)).^2;
A=RaisedCosineNormSquare(1:(length(RaisedCosineNormSquare)/2));
rampUp = A;
rampDown = fliplr(rampUp);
midT1 = ones(1,length(t1) - length(rampUp) - length(rampDown));
envelopeT1 = [rampUp midT1 rampDown];
midT2 = ones(1,length(t2) - length(rampUp) - length(rampDown));
envelopeT2 = [rampUp midT2 rampDown];
pad = zeros(1,50);

Vol1 = t1.*envelopeT1;
Vol2 = t2.*envelopeT2;

temp1 = [Vol1 zeros(size(t_beep))]; % Create adapter & tester tones.
temp2 = [Vol2 zeros(size(t_beep))];

tones(1,:) = [repmat(temp1,1,11) repmat(zeros(size(t_beep)),1,8)];
tones(2,:) = [repmat(temp1,1,11) repmat(zeros(size(t_beep)),1,2) repmat(temp1,1,3)];
tones(3,:) = [repmat(temp1,1,11) repmat(zeros(size(t_beep)),1,2) repmat(temp2,1,3)];

screenNumber = max(Screen('Screens')); %get screen number
AssertOpenGL; %Initalise

% Get the size of the display (rect) and the id of the new window (w).
[w, rect] = Screen('OpenWindow', screenNumber, [140 140 140]);

% Enable alpha blending with proper blend-function. We need it
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[center(1), center(2)] = RectCenter(rect);
fps=Screen('FrameRate',w);      % frames per second
ifi=Screen('GetFlipInterval', w);
if fps==0
   fps=1/ifi;
end;

black = [140 140 140]; %grey
white = [0 161 0]; %green
HideCursor;	% Hide the mouse cursor

Priority(MaxPriority(w)); %priortize the window we want

% Set text size:
Screen('TextSize', w, 30);
Screen('TextStyle',w,1);

% Do initial flip...
Screen('Flip', w);

% Calculate general parameters for display:
ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
fix_cord = [center-fix_r*ppd center+fix_r*ppd];         % coords for fixation point (center)
countr = 0;
ttt = 1;
for runs = 1:3
    
        % Randomize the adaptation order:
        adapt_order = Shuffle(1*ones(num_trials,1)); %FOR CURRENT PURPOSES: all stimuli are the same (hence all ones)
        
        % Randomize the iti:
        iti_num = ceil(num_trials/3);
        ii = 1;
        for n_iti = 1:num_trials
            iti_value(ii,1) = (rand*4+6);
            ii = ii+1;
        end
        iti_order = iti_value.*ones(num_trials,1);
        
        % Number of letters to show at fixation:
        num_letters = (1/letter_length)*sum(iti_order)+100+3.15*num_trials+2000;
        
        % Shuffle their order (KbName 65 to 90 is entire alphabet):
        letter_order = Shuffle(repmat([76:90],1,200));
        letter_order = letter_order(1:num_letters);
        for nn = 1:num_letters
            if nn > length(letter_order)
                break
            end
            templetter = letter_order(nn);
            if nn+1<=num_letters
                if nn+1 > length(letter_order)
                    break
                end
                if templetter == letter_order(nn+1)
                    if nn+13 <=num_letters
                        for nnn = nn+1:nn+11
                            if nnn+1 > length(letter_order)
                                break
                            end
                            if letter_order(nnn) == letter_order(nnn+1)
                                letter_order(nnn+1) = [];
                            end
                        end
                    end
                end
            end
        end

    if eyetracking == 1
            dummymode = 0;
            el = EyelinkInitDefaults(w);
            el.backgroundcolour = [140 140 140];
            el.foregroundcolour = [0 161 0];
            el.calibrationtargetcolour = [0 161 0];
            EyelinkUpdateDefaults(el);
            if ~EyelinkInit(dummymode)
                fprintf('Eyelink Init aborted.\n');
                Eyelink('Shutdown');
                sca;
                commandwindow;
                ListenChar(0);
                return;
            end

            [v vs] = Eyelink('GetTrackerVersion');
            fprintf('Running experiment on a ''%s'' tracker.\n',vs);
            if sense == 0
                if runs == 1
                    edfFile = sprintf(strcat(num2str(number),'.edf'));
                elseif runs == 2
                    edfFile = sprintf(strcat(num2str(number),'v2.edf'));
                else
                    edfFile = sprintf(strcat(num2str(number),'v3.edf'));
                end
            else
                if runs == 1
                    edfFile = sprintf(strcat(num2str(number),'a.edf'));
                elseif runs == 2
                    edfFile = sprintf(strcat(num2str(number),'av2.edf'));
                else
                    edfFile = sprintf(strcat(num2str(number),'av3.edf'));
                end
            end
            Eyelink('Openfile',edfFile);

            EyelinkDoTrackerSetup(el);
            EyelinkDoDriftCorrection(el);
            Eyelink('StartRecording');
            WaitSecs(0.5);

    end

    %%% Run the experiment:

    buttonKeys = [32]; % Set button keys (response box)
    escapeKey = KbName('ESCAPE'); % Set escape key
    ListenChar(2); % Prevent pressed keys from showing up in command prompt
    KbWait([], 2);
    dur0 = GetSecs;
    terminate = 0;  % Change to 1 when stopping experiment.
    if eyetracking == 1
        Screen('DrawText',w,'Please look at the circle that will appear on the',30,250,white);
        Screen('DrawText',w,'screen for the next minute, until it disappears.',30,300,white);
        Screen('Flip',w);
        pause(10)
        Screen('FillOval',w,uint8([0 161 0]),fix_cord);
        Screen('Flip',w);
        Eyelink('Message','Baseline Recording Start');
        WaitSecs(45);
        Eyelink('Message','Baseline Recording End');
    end
    Screen('DrawText',w,'In the following game, press the space bar when you',30,250,white);
    Screen('DrawText',w,'see a letter appear twice in a row. Please sit still',30,300,white);
    Screen('DrawText',w,'and keep your eyes at the center of the screen.',30,350,white);
    Screen('DrawText',w,'Press the right arrow key to begin.',30,450,white);
    Screen('Flip',w);
    KbWait([],2);

    % Start with n seconds of fixation only:
    for i = 1:15*fps

        % Show fixation:
        Screen('FillOval', w, uint8([0 161 0]), fix_cord);
        Screen('Flip', w);

        % Terminate if escape key is pressed
        [keyIsDown, seconds, keyCode] = KbCheck;

    end
    frame_count = 1; % Use to count frames between letters
    letter_count = 1; % Use to keep track of which letter to display
    response_count = 1; % Use to keep track of responses
    response_correct = 0; % Use to keep track of button press feedback
    oval_count = 0; % Keep track of how long oval color was displayed
    responses = []; % Create empty responses matrix
    for trial = 1:num_trials

        if terminate == 1   % If escape key was pressed, break loop.
            break;
        end

        %%%% Stimulus (adaptor & tester) and iti %%%%%%

        % Send message to eye tracker that trial is starting (adapter&test):

        if eyetracking == 1
            Eyelink('Message','Adapter is starting');
            if adapt_order(trial) == 1
                Eyelink('Message','No test');
            elseif adapt_order(trial) == 2
                Eyelink('Message','Adapter trial');
            else
                Eyelink('Message','Test trial');
            end
        end
        t_start = GetSecs;
        s = [tones(adapt_order(trial),:); tones(adapt_order(trial),:)]';
        sound(s,freq);

        % Iter trial interval:
        while GetSecs <= t_start+iti_order(trial)+3.15

            % Auditory stimulus ends after 3.15 seconds. Tell eye tracker
            % that stimulus has ended at that point:
             if GetSecs == 3.15+t_start
                if eyetracking == 1
                    Eyelink('Message','Adapter has ended.');
                end
             end
            % Draw the fixation:
            if response_correct == 1
                oval_count = oval_count + 1;
            elseif response_correct == 2
                oval_count = oval_count + 1;
            else
                oval_count = 0;
            end
            Screen('FillOval', w, uint8(white), fix_cord);
            if oval_count == 0.5*fps    % Stop displaying feedback after half a second
                response_correct = 0;
                oval_count = 0;
            end
            % Draw a letter at fixation:
            if frame_count < letter_length*fps-5    % Leave a gap between letters of 5 frames
                if response_correct == 1
                    Screen('DrawText', w, KbName(letter_order(letter_count)), center(1)-13, center(2)-25, [124 124 255]);
                elseif response_correct == 2
                    Screen('DrawText',w,KbName(letter_order(letter_count)),center(1)-13,center(2)-25,[255 63 0]);
                else
                    Screen('DrawText', w, KbName(letter_order(letter_count)), center(1)-13, center(2)-25, black);
                end
            end
            % Advance to the next letter when reaching letter length
            if frame_count == letter_length*fps
                frame_count = 1;
                letter_count = letter_count+1;
            else
                frame_count = frame_count + 1;  % Advance frame_count
            end
            % When the same letter appears repeatedly (one back):
            if frame_count == 1 && letter_order(letter_count) == letter_order(letter_count-1)
                start_time = GetSecs;
                if eyetracking == 1
                Eyelink('Message','Letter repeated');
                end
                countr = countr+1;
                if ttt == 1
                    timere(ttt) = GetSecs;
                    ttt = 2;
                else
                    timere(ttt) = GetSecs-sum(timere(1:ttt-1));
                    ttt = ttt+1;
                end
            end
            % When a key is pressed:
            [keyIsDown, press_time, keyCode, deltaSecs] = KbCheck;
            if keyIsDown
                if keyCode(escapeKey)
                    terminate = 1;
                    break;
                elseif sum(keyCode(buttonKeys))>0
                    if exist('start_time') && press_time-start_time < response_time
                        if oval_count < 1
                            response_correct = 1;   % correct button press
                            responses(response_count,1) = press_time-start_time;
                            responses(response_count,2) = double(find(keyCode(buttonKeys)==1));
                            responses(response_count,3) = letter_count;
                            oval_count = 1;
                            if eyetracking == 1
                            Eyelink('Message','Correct');
                            end
                        end

                    else
                        if oval_count < 1
                            response_correct = 2;   % incorrect button press or too late
                            responses(response_count,1) = 0;
                            responses(response_count,2) = double(find(keyCode(buttonKeys)==1));
                            responses(response_count,3) = letter_count;
                            oval_count = 1;
                            if eyetracking == 1
                            Eyelink('Message','Wrong');
                            end
                        end
                    end
                    response_count = response_count + 1;
                    clearvars keyIsDown press_time keyCode deltaSecs KbCheck
                end
            end
            Screen('Flip', w);

        end

    end

    dur = GetSecs-dur0;
    extra = ceil(315-dur);

    % Save experiment structure and responses:
    tones_hz = [freq_t1 freq_t2];
    if runs == 1
        save([path sub num2str(number) '_Order_Aud_' datestr(now, 'yymmdd_HHMMSS')], 'adapt_order', 'iti_order', 'tones_hz', 'letter_order', 'dur');
        save([path sub num2str(number) '_Responses_Aud_' datestr(now, 'yymmdd_HHMMSS')], 'responses');
    else
        save([path sub num2str(number) 'v2_Order_Aud_' datestr(now,'yymmdd_HHMMSS')],'adapt_order','iti_order','tones_hz','letter_order','dur');
        save([path sub num2str(number) 'v2_Responses_Aud_' datestr(now,'yymmdd_HHMMSS')],'responses');
    end
    % End with n seconds of fixation only:

    for i = 1:extra*fps
        if terminate == 1   % If escape key was pressed, break loop.
            break;
        end        
        % Show fixation spot:
        Screen('FillOval', w, uint8(white), fix_cord);	% draw fixation dot (flip erases it)
        Screen('Flip', w);

        % Terminate if escape key is pressed:
        [keyIsDown, seconds, keyCode] = KbCheck;
        if keyCode(escapeKey)
            break;
        end
    end

    if eyetracking == 1
            Eyelink('StopRecording');
            Eyelink('CloseFile');
            Eyelink('Shutdown');
    end
    if runs < 3
        Screen('DrawText',w,'You may now get up and take a break.',30,250,white);
        Screen('DrawText',w,'Press the right arrow key to continue.',30,300,white);
        Screen('Flip',w);

        KbWait([],2); 
    end
        
end
    
sca;
commandwindow;
ListenChar(0);

Priority(0);
totdur = GetSecs-dur0;
ShowCursor
Screen('CloseAll');    
ListenChar(0);

Priority(0);
ShowCursor
Screen('CloseAll');
ListenChar(0);  % Return writing capabilities   
