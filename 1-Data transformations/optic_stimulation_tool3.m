function optic_stimulation_tool2()
% OPTIC_STIMULATION_TOOL2()
% A command line script to control the optical stimulation.
% This is the 2nd version of this stimulator. It is tailored to mix and compare different kinds of stimuli.
%
% Commands:
%   h       this help
%   end     exit script
%
%   t 1/0    switch test mode on/off
%   w 1/0    1 for white on black; 0 for black on white
%   q 1/0    switch between "queue" mode and "24 cycle" mode           
%
%   TYPES OF STIMULI:
%       0    Empty (nothing)
%       f    full flash
%       l    flashing
%       c/C  crash (1s) and fast crash (0.5 s)
%       o/O  slow crash (1.5 and 2 s)
%       u    superslow crash (3s)
%       i/I  robotic ringy crash, normal speed (1s) and fast (0.5 s)
%       y/Y  robotic ringy crash, slow (0.5 and 2 s)
%       1-5  robotic rings (from the robotic crash) added one by one. Aren't removed from the screen!
%       b    backwards crash (a ring contracting from the sides to the center)
%       r    ramp
%       g    grid (alternative to crash and ramp)
%       p    psycho (expading rings)
%       s    stripes changing color, similar to Ed. Ruthazer's experiments
%       !    change background color
%       a number - adjust the x-position of the stimulus (in case the rig misbehaves). NOT FULLY IMPLEMENTED

% Apr 25 2012: Rev.
% Jun 13 2013: Adjusting x position.
% Jan 07 2014: + Digitized (downsampled) crash, and reshuffled (scrambled) crash.
%               At the same time, some weird "crashes" (C, O) are removed, epileptogenic_grid renamed to 'u',
%               and 's' now stands for "Scrambled".

fprintf('\n]]] === OPTIC_STIMULATION_TOOL started\n');
clear;

%%% ========================== CONSTANTS =============================

%%% -------- Mode holders (may be changed later in the cycle
test_mode = 1;          % In test mode it does not wait for Clampex, and it does not make wasted sweeps
wait_in_test_mode = 3;  % How many seconds to wait between stimulus presentations in test mode
whiteOnBlack = 0;       % 1 for classic white-on-black, 0 for reverse black-on-white
flagMainCycle = 1;      % when this is zeroed, the program would stop
flagQueue = 1;          % If 1 - follows the program once; if 0 - fills 24 repetitions with this program

%%% ------- Universal constants
whichScreen = 2;    % 0 for default; 1 for main (which is the same); our choice is 2
screenResolution = [1024 768];
screenWidth = 150;  % Screen Size in pixels. 100 square fully fits in the fiber, 120-150 square fully covers it
nCircles = 5;
xgrid = 5;
xAdj = 0;
nPixels = 17;       % For downsampled stimuli. We agreed that 17 is a good number

stimulusDuration = 5;       % Total duration. 5 for physiology, 2 for testing
dynamicDurationBase = 1;    % How quickly crashes and ramps grow - may be overriden later
flashRate = 0.25;           % Flashing phase (either on or off, not total) duration in ms
nFrames = 50;
flagHide = 1;               % Whether the screen should be cleared after the stimulus presentation
nFlashes = floor(stimulusDuration/flashRate);
[x,y] = meshgrid(1:xgrid);
squareEdges = [x(:)-1 y(:)-1 x(:)-1 y(:)-1]'*screenWidth/xgrid; % Useful alias
v = screenWidth/2;          % Just a useful alias
myrect = repmat(screenResolution/2,1,2) + screenWidth/2*[-1 -1 1 1];
                % rect: left top right bottom.
                
mySound = sin((1:1000)/10000*2*pi*440*4); % Short beep
                
%%% -- This part is about reshuffling
lowR_m = (nPixels+1)/2;                   % Middle index for the downsampled matrix
[lowR_x,lowR_y] = meshgrid(1:nPixels);     % Grid of coordinates
crashMatrix = sqrt((lowR_x-lowR_m).^2+(lowR_y-lowR_m).^2);  % Distances
crashMatrix = crashMatrix/crashMatrix(1,1);                         % Normalized distances
scrmbMatrix = crashMatrix;                          % Future reshuffled matrix
scrmbMatrix(:) = crashMatrix(randperm(nPixels^2));     % Actual act of reshuffling


%%% ========================== GLOBAL INITIALIZATION =============================
rand('seed',floor(sum(clock)*1000)); % Random generator - otherwise all the experiments get the same
filename = [datestr(now(),'yyyymmdd') '-StimLog.txt'];
fH = fopen(filename,'a');
fprintf(fH,['\r\nV2 Initialized: ' datestr(now()) '\r\n']);

fprintf('   Suppressing warnings from both Screen and IOPort commands.\n');
screenDebugLevel = Screen('Preference','VisualDebugLevel',0);
screenVerbosityLevel = Screen('Preference', 'Verbosity',0);
ioPortVerbosityLevel = IOPort('Verbosity',0);
Screen('Preference', 'SkipSyncTests', 1); fprintf('   Bypassing SyncTest.\n');
%     % The Psychbox recommends the 'SkipSyncTests' instruction as a 
%     % risky measure to bypass the synchronization test,
%     % which sometimes stopped this script, returning an error.
window = Screen('OpenWindow', whichScreen);

% Usage: [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = 
%   Screen('Flip', windowPtr [, when] [, dontclear] [, dontsync] [, multiflip]);
% The interesting part here is "dontsync" which is =0 by default, but can be set to 2 (draw immediately)
% Timestamps mean actually: beginning, stimulus onset, out of flip. 3-1=flip duration.
Screen(window, 'FillRect', 255*(1-whiteOnBlack)); % Make the screen of background color
Screen(window, 'Flip');
if(test_mode)
    fprintf('WORKING in TEST MODE. Use t command to change this.\n');
end
fprintf('   You can enter commands now.\n');

%%% ========================== MAIN CYCLE =============================
while(flagMainCycle)
    flagStartExperiment = 0; % Will be set to 1 as soon a protocol is initiated
    while(~flagStartExperiment)
        %%% --------------------- Get input and analyze it ----------------------
        command = input(']] ','s');        
        if(sum(command==' ')>0) % Commands with modifier
            temp = find(command==' ');
            modifier = command(temp(1)+1:end);
            command = command(1:temp(1)-1); % Cut the tail
        else
            modifier = '';
        end
        
        switch command
            case 'end'
                flagMainCycle = 0;
                command = ''; % A 'secret command' to skip experimental part
                break;
            case {'help','h'}
                help optic_stimulation_tool2;
            case {'test','t'}
                if(~isempty(modifier))
                    if(ismember(str2num(modifier),[0 1]))
                        test_mode = str2num(modifier);
                        fprintf('Test mode set to %d\n',test_mode);
                    else
                        fprintf('Test mode constant not recognized\n');
                    end
                else
                    fprintf('Test mode = %d\n',test_mode);
                end
            case 'a'
                if(isempty(modifier))
                    fprintf('Current adjustment = %d\n',xAdj);
                end
            case 'w'
                if(~isempty(modifier))
                    if(ismember(str2num(modifier),[0 1]))
                        whiteOnBlack = str2num(modifier);
                        fprintf('WhiteOnBlack mode flag set to %d\n',whiteOnBlack);
                    else
                        fprintf('WhiteOnBlack mode flag not recognized\n');
                    end
                    Screen(window, 'FillRect', 255*(1-whiteOnBlack)); % Make the screen black
                    Screen(window, 'Flip');    
                else
                    fprintf('WhiteOnBlack mode flag = %d\n',whiteOnBlack);
                end
            case {'queue','q'}
                if(~isempty(modifier))
                    if(ismember(str2num(modifier),[0 1]))
                        flagQueue = str2num(modifier);
                        fprintf('Queue mode is set to %d\n',flagQueue);
                    else
                        fprintf('Queue mode not recognized\n');
                    end
                else
                    fprintf('Queue mode = %d\n',flagQueue);
                end
            case ''
                % Nothing
            otherwise
                % May be a protocol formula
                goodProtocol = 1;
                for(iChar=1:length(command))
                    switch(command(iChar))
                        case {'f','l','c','b','r','g','p','s','0','!','o','i','y',...
                                'u','1','2','3','4','5'}
                            % Fine
                        otherwise
                            goodProtocol = 0;
                            fprintf('What the heck would %s mean???\n',command(iChar));
                    end
                end
                if(goodProtocol) % All the caracters proved to be legal
                    flagStartExperiment = 1;
                else
                    fprintf('Command is not recognized \n');
                end
        end
    end
        
    %%% =============== LOCAL INITIATION ============    
    if(isempty(command)) % Termination sequence
        break;
    end
    
    if(flagQueue)
        nRepetitions = length(command);        
    else
        nRepetitions = 12;
        fprintf('Using %d repetitions, mind it!\n',nRepetitions);
    end
    
    % Setting the colors. Below I just assign 0 and 255 values as black and white respectively,
    % but principally a more elaborate way of doing this is through WhiteIndex(window) and
    % BlackIndex functions. It may make the program less platform-specific I guess...
    if(whiteOnBlack)
        colorFor = 255;
        colorBak = 0;
    else % Inverted picture
        colorFor = 0;
        colorBak = 255;
    end
    
    % To be secured agains trains of commands I open and close IOPort each time
    [handle_com, errmsg] = IOPort('OpenSerialPort', 'COM1', 'ReceiveTimeout=100.0');
    IOPort('Flush', handle_com);

    %%% ============== REPORTING =================
    fprintf(fH,['\r\n Protocol: %s . ' datestr(now()) '\r\n'],command);
    fprintf(fH,'OnBlack %d, WidthPX %d, Grid %d, Repts %d, StimDur %d \r\n',...
        whiteOnBlack,screenWidth,xgrid,nRepetitions,stimulusDuration);
    
    fprintf('   Protocol: %s. ',command);
    fprintf('OnBlack %d, WidthPX %d, Grid %d, Repts %d, StimDur %d \r\n',...
        whiteOnBlack,screenWidth,xgrid,nRepetitions,stimulusDuration);
    
    fprintf(']]] Start the Clampex now!\n');
    ding;

    Screen(window, 'FillRect', colorBak); % Set the background
    Screen(window, 'Flip'); % Show the black screen
    ticCycle = tic; % Just for the toc below not to produce error    
    
    for(iStim=1:nRepetitions) % -------------------- Repetitions cycle        
        experimentType = command(mod(iStim-1,length(command))+1); % Rotate through command
        flagGoFurther = waitForNextSweep(test_mode,wait_in_test_mode,handle_com);
        if(~flagGoFurther)
            break
        end        
        ticCycle = tic; % Just for the toc below not to produce error
        fprintf('Stimulus #%2d/#%2d: %s - on',iStim,nRepetitions,experimentType);
        dynamicDuration = dynamicDurationBase; % Default value
        frameDelay = dynamicDuration/nFrames; % Also a default value - needs to be recalculated if dD was changed
        
        %%% ========================== MAIN SWITCH =============================
        switch experimentType    
            case '0'
                %%% ============ EMPTY SLOT ================
                WaitSecs(stimulusDuration); % Show the stimulus for that long
                % WaitSecs is the same like PAUSE, only from the PsychToolbox, and so presumably more efficient
                
            case '!'
                %%% ============ CHANGE BACKGROUND COLOR ================
                whiteOnBlack = 1-whiteOnBlack;
                Screen(window, 'FillRect', 255*(1-whiteOnBlack)); % Make the screen background
                Screen('Flip',window,0,0,2); % Show the picture
                WaitSecs(stimulusDuration); % Show the stimulus for that long
                if(whiteOnBlack)
                    colorFor = 255;     colorBak = 0;
                else
                    colorFor = 0;       colorBak = 255;
                end
                
            case 'f' %%% ========================== FULL FIELD FLASH =============================            
                Screen(window, 'PutImage', colorFor*ones(screenWidth));
                Screen('Flip',window,0,0,2); % Show the picture
                WaitSecs(stimulusDuration); % Show the stimulus for that long

            case 'l' %%% ========================== EPILEPTOGENIC FULL FIELD FLASHING =============================
                currentPhase = 1;
                
                for(q=1:nFlashes)
                    if(currentPhase)
                        Screen('FillRect', window, colorFor, myrect);
                        fprintf('1');
                    else
                        Screen('FillRect', window, colorBak, myrect);
                        fprintf('0');
                    end
                    currentPhase = 1-currentPhase;
                    Screen('Flip',window,0,0,2); % Show the picture
                    WaitSecs(flashRate); % Show the stimulus for that long                    
                end
                                
            case 'p' %%% ========================== EPILEPTOGENIC CIRCLES =============================
                while(toc(ticCycle)<=stimulusDuration)
                    phase = toc(ticCycle)/(flashRate*2); % Continous
                    phase = phase-floor(phase);        % 0 to 1 saw-shape
                    for(q=1:nCircles)
                        r = screenWidth/2/nCircles*(nCircles-q+phase*2);
                        colorThis = mod(q,2)*255;  % WhiteOnBlack doesn't matter for this stimulus
                        Screen('FillOval', window, colorThis, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*r, screenWidth);
                                % rect: left top right bottom.                                                
                    end
                    Screen('Flip',window,0,0,2); % Show the picture
                end                               

            case 'u' %%% ========================== EPILEPTOGENIC GRID =============================                
                while(toc(ticCycle)<=stimulusDuration)                    
                    phase = tri(toc(ticCycle)/flashRate); % Triangle wave
                    for(q=1:nCircles)                        
                        r = screenWidth/2/nCircles*(nCircles-q+1);
                        colorThis = 255*(mod(q,2)*phase+(1-mod(q,2))*(1-phase)); % WhiteOnBlack doesn't matter for this stimulus
                        Screen('FillOval', window, colorThis, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*r, screenWidth);
                                % rect: left top right bottom.                                                
                    end
                    Screen('Flip',window,0,0,2); % Show the picture
                end

            case {'c','o','s'} %%% ========================== CRASH (Low-res) =============================
                switch experimentType
                    case 'c'; dynamicDuration = 1; scrambleFlag = 0;
                    case 'o'; dynamicDuration = 2; scrambleFlag = 0;
                    case 's'; dynamicDuration = 1; scrambleFlag = 1;
                end
                frameDelay = dynamicDuration/nFrames;
                ticStim = tic;
                while(toc(ticStim)<=stimulusDuration)                    
                    if(~scrambleFlag)
                        % normal crash
                        Screen('PutImage', window, imresize(255*(crashMatrix>=(toc(ticStim)/dynamicDuration)),screenWidth/nPixels,'nearest'));
                    else
                        % scrambled crash
                        Screen('PutImage', window, imresize(255*(scrmbMatrix>=(toc(ticStim)/dynamicDuration)),screenWidth/nPixels,'nearest'));
                        % scrmbMatrix
                    end
%                     Screen('FillOval', window, colorFor, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*radius, screenWidth);
%                         % rect: left top right bottom.                                                                    
                    Screen(window, 'Flip'); % Now show the picture                    
                end                
                
%             case {'nope'} %%% ========================== OLDschool REAL CRASH =============================
%                 switch experimentType
%                     case 'c'; dynamicDuration = 1;
%                     case 'C'; dynamicDuration = 0.5;
%                     case 'o'; dynamicDuration = 1.5;
%                     case 'O'; dynamicDuration = 2;
%                     case 'u'; dynamicDuration = 3;
%                 end
%                 frameDelay = dynamicDuration/nFrames;                
%                 for(q=1:nFrames)                                    
%                     frameTic = tic;
%                     radius = v*q/nFrames;
%                     Screen('FillOval', window, colorFor, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*radius, screenWidth);
%                         % rect: left top right bottom.                                                                    
%                     Screen(window, 'Flip'); % Now show the picture                            
%                     WaitSecs(frameDelay-toc(frameTic)); % Show the frame for that long                    
%                 end
%                 WaitSecs(stimulusDuration-dynamicDuration);
                                
            case 'b' %%% ========================== BACK CRASH =============================                            
                for(q=1:nFrames)                
                    frameTic = tic;
                    Screen('FillOval', window, colorFor, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*v, screenWidth);
                        % Full screen
                    radius = v*sqrt(1-(q/nFrames)^2); % And this one will contract
                    Screen('FillOval', window, colorBak, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*radius, screenWidth);
                        % rect: left top right bottom.                                                                    
                    Screen(window, 'Flip'); % Now show the picture                            
                    WaitSecs(frameDelay-toc(frameTic)); % Show the frame for that long
                end
                WaitSecs(stimulusDuration-dynamicDuration);
            
            case {'i','I','y','Y'} %%% ========================== ROBOTIC CRASH =============================            
                localNFrames = 5;
                switch experimentType
                    case 'i'; dynamicDuration = 1;
                    case 'I'; dynamicDuration = 0.5;
                    case 'y'; dynamicDuration = 1.5;
                    case 'Y'; dynamicDuration = 2;
                end
                frameDelay = dynamicDuration/localNFrames;                                
                stimStartTic = tic;
                for(q=1:localNFrames)
                    frameTic = tic;
                    radius = v*q/localNFrames;
                    Screen('FillOval', window, colorFor, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*radius, screenWidth);
                        % rect: left top right bottom.                                                                    
                    Screen(window, 'Flip'); % Now show the picture                            
                    WaitSecs(frameDelay-toc(frameTic)); % Show the frame for that long                    
                end
                WaitSecs(stimulusDuration-dynamicDuration);
                
            case {'1','2','3','4','5'} %%% ========================== ROBOTIC CRASH DISASSEMBLED =============================            
                localNFrames = 5;
                ringNumber = eval(experimentType);                
                radius = v*ringNumber/localNFrames;
                Screen('FillOval', window, colorFor, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*radius, screenWidth);
                    % rect: left top right bottom.                                                                    
                Screen(window, 'Flip'); % Now show the picture                                            
                WaitSecs(stimulusDuration-dynamicDuration);
                if(ringNumber<5)
                    flagHide = 0; % Don't clear this stimulus off
                end
            
            case 'r' %%% ========================== COLOR RAMP =============================            
                for(q=1:nFrames)                
                    frameTic = tic;
                    if(whiteOnBlack)
                        currentColor = q^2/nFrames^2*255;
                    else
                        currentColor = (1-q^2/nFrames^2)*255;
                    end
                    Screen('FillOval', window, currentColor, repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*v, screenWidth);
                    Screen(window, 'Flip'); % Now show the picture                            
                    WaitSecs(frameDelay-toc(frameTic)); % Show the frame for that long
                end
                WaitSecs(stimulusDuration-dynamicDuration);

            case 'g' %%% ========================== GRID (FORMER FLOWER) =============================
                for(q=1:nFrames)                
                    frameTic = tic;                                                            
                    Screen('FillRect', window, colorFor,...
                            repmat(repmat(screenResolution,1,2)'/2+[-1 -1 1 1]'*v/xgrid*q/nFrames,1,xgrid^2)+squareEdges...
                            + v/xgrid/2 - v,...
                            screenWidth);
                        % To fill n rectangles, provide "rect" as a 4 rows by n columns matrix                    
                    Screen(window, 'Flip'); % Now show the picture                            
                    WaitSecs(frameDelay-toc(frameTic)); % Show the frame for that long
                end
                WaitSecs(stimulusDuration-dynamicDuration);
                
            case 'end'
                % Nothing
            otherwise
                fprintf('Possible bug: Protocol was not found. \n');
        end % ----------- Experiment type switch
        
        if(flagHide)
            Screen(window, 'FillRect', colorBak); % Set the background back to background =)
            Screen(window, 'Flip');
            fprintf(' - off.\n');        
        else
            fprintf(' (kept)\n');
        end
        flagHide = 1;
        %sound(mySound,10000); % Make some sound        
        if((nRepetitions==24)&&(iStim==12))
            ding('beep2');
        else
            ding('beep');
        end
        
    end % ----------- Repetitions cycle
    
    % == Local wrap-up (entering standby mode) ==
    IOPort('CloseAll');    
    fprintf(fH,['Protocol %s completed: ' datestr(now()) '\r\n'],command);    
    ding;
end % ----------- Main cycle

%%% ========================== GLOBAL WRAP-UP =============================
WaitSecs(1); % For all possible Clampex recordings to finish
Screen('CloseAll');

Screen('Preference','VisualDebugLevel',screenDebugLevel);
Screen('Preference', 'Verbosity',screenVerbosityLevel);
Screen('Preference','VisualDebugLevel',screenVerbosityLevel);

fprintf(fH,['Terminated: ' datestr(now()) '\r\n\r\n']);
fclose(fH);
fprintf(']]] === Program terminated\n');

end

function res = waitForNextSweep(test_mode,wait_in_test_mode,handle_com);
% Either wait for Clampes, or (in test mode) - just pause a bit
res = 1;
[~,~,key] = KbCheck;
escapeKey = KbName('esc');
if(key(escapeKey))
    res = 0;
    return
else
    if(test_mode)
        WaitSecs(wait_in_test_mode); fprintf(' Test mode: pausing for %2.1f sec. ',wait_in_test_mode);    
    else
        IOPort('Read', handle_com, 1 , 1); % Wait for the trigger
    end
end
end

function b = tri(a)
% Triangle function with period of 2
a = mod(a,2);
b = a.*(a<=1)+(2-a).*(a>1);
end