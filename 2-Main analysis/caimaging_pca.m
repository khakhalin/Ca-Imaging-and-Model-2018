function caimaging_pca()
% caimaging_pca()
%
% Analysis of network responses.

% Oct 05 2017: Caimaging_pca forked from caimaging_structure. All TE calculations and other heavy lifting is left there,
%   while this program now does everything else. Basically, all analysis that doesn't use connectivity reconstruction.
% Aug 14 2018: Lots of minor improvements.
% Oct 01 2018: General clean-up

S = [];
iBrain = 0;

%%% Base input folders are specified below as they are different for s46 and s49 experiments
auxFolder = 'C:\_Data\___Ca imaging\auxData\';   % A folder to save some selected results. It's the same folder that contains reconstructions.

% % %%% --- stage 46 set
iBrain = iBrain+1; folderName{iBrain} = '140718a'; age(iBrain) = 46;
iBrain = iBrain+1; folderName{iBrain} = '140716b'; age(iBrain) = 46;
iBrain = iBrain+1; folderName{iBrain} = '140716a'; age(iBrain) = 46;
iBrain = iBrain+1; folderName{iBrain} = '140715'; age(iBrain) = 46;
iBrain = iBrain+1; folderName{iBrain} = '140711'; age(iBrain) = 46;

iBrain = iBrain+1; folderName{iBrain} = '140708b'; age(iBrain) = 46;
iBrain = iBrain+1; folderName{iBrain} = '140708a'; age(iBrain) = 46;
iBrain = iBrain+1; folderName{iBrain} = '140705a'; age(iBrain) = 46;
iBrain = iBrain+1; folderName{iBrain} = '140704b'; age(iBrain) = 46; % strange responses, weird PCA. Two cells are very weird. Check.
iBrain = iBrain+1; folderName{iBrain} = '140704a'; age(iBrain) = 46; % Weak
iBrain = iBrain+1; folderName{iBrain} = '140627'; age(iBrain) = 46;
iBrain = iBrain+1; folderName{iBrain} = '140626';  age(iBrain) = 46;
%%% iBrain = iBrain+1; folderName{iBrain} = '140620'; age(iBrain) = 46; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140619b'; age(iBrain) = 46; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140619a'; age(iBrain) = 46; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140613'; age(iBrain) = 46; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140611'; age(iBrain) = 46; % excluded
iBrain = iBrain+1; folderName{iBrain} = '140610'; age(iBrain) = 46;
%%% iBrain = iBrain+1; folderName{iBrain} = '140530b'; age(iBrain) = 46; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140530a'; age(iBrain) = 46; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140529'; age(iBrain) = 46; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140528'; age(iBrain) = 46; % excluded
iBrain = iBrain+1; folderName{iBrain} = '140516'; age(iBrain) = 46;
%%% iBrain = iBrain+1; folderName{iBrain} = '140502'; age(iBrain) = 46; % excluded

% % %%% --- stage 49 set
iBrain = iBrain+1; folderName{iBrain} = '140726'; age(iBrain) = 49;
iBrain = iBrain+1; folderName{iBrain} = '140724'; age(iBrain) = 49; % pretty retinotopy
iBrain = iBrain+1; folderName{iBrain} = '140723'; age(iBrain) = 49;
iBrain = iBrain+1; folderName{iBrain} = '140722'; age(iBrain) = 49; % Was used for raw data images
iBrain = iBrain+1; folderName{iBrain} = '140718b'; age(iBrain) = 49; % pretty retinotopy
iBrain = iBrain+1; folderName{iBrain} = '140714'; age(iBrain) = 49; % Late responses (rebounds) for crash
iBrain = iBrain+1; folderName{iBrain} = '140710'; age(iBrain) = 49;
iBrain = iBrain+1; folderName{iBrain} = '140709'; age(iBrain) = 49; % pretty retinotopy - the best
iBrain = iBrain+1; folderName{iBrain} = '140707'; age(iBrain) = 49;
iBrain = iBrain+1; folderName{iBrain} = '140705b'; age(iBrain) = 49;
%%% iBrain = iBrain+1; folderName{iBrain} = '140703'; age(iBrain) = 49; % excluded
iBrain = iBrain+1; folderName{iBrain} = '140612'; age(iBrain) = 49;
%%% iBrain = iBrain+1; folderName{iBrain} = '140522'; age(iBrain) = 49; % excluded
iBrain = iBrain+1; folderName{iBrain} = '140521'; age(iBrain) = 49;
%%% iBrain = iBrain+1; folderName{iBrain} = '140505'; age(iBrain) = 49; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140408'; age(iBrain) = 49; % excluded
iBrain = iBrain+1; folderName{iBrain} = '140328'; age(iBrain) = 49; % weak, bad retinotopy
%%% iBrain = iBrain+1; folderName{iBrain} = '140326'; age(iBrain) = 49; % excluded, cfs2 is shorter than others
%%% iBrain = iBrain+1; folderName{iBrain} = '140325'; age(iBrain) = 49; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140318'; age(iBrain) = 49; % excluded
%%% iBrain = iBrain+1; folderName{iBrain} = '140317'; age(iBrain) = 49; % excluded
iBrain = iBrain+1; folderName{iBrain} = '140314'; age(iBrain) = 49; % weak, bad retinotopy
%%% iBrain = iBrain+1; folderName{iBrain} = '140312'; age(iBrain) = 49; % ex cluded
iBrain = iBrain+1; folderName{iBrain} = '140311'; age(iBrain) = 49; % Noisy, except for 2 cells, but produced a decent graph
iBrain = iBrain+1; folderName{iBrain} = '140310'; age(iBrain) = 49; % Short recording


% ------------- Control settings and constants -------------
showRawData = 0;                    % Whether raw (well, raw-ish) data debugging figures need to be shown
showRawFigure = 0;                  % Not all raw data, but just enough for a figure (averages for this brain and a band of CI around each average)
pcaStimType = 1;                    % Set to 1 if PCA and response dynamics are to be run on crashes only (recommended); 2=flash, 3=scramble; 0=full triad
randomizePositions = 0;             % Set to 1 if all xy positions should be randomly reassigned (as a H0-style test for position-related methods)
showAverageResponses = 0;           % Show average responses
flagUsePeakAmps = 0;                % If set to 1, overrides cumulative amplitudes (amps) with peak amplitudes (ampsPeak). NOTE that it rewrites selFC in AUXDATA
peakRange = [-0.1 0.1];             % Range around the average flash peak at which amplitudes should be taken

doPCA = 1;                          % Factor analysis. The "origin point" for retinotopic responses is also calculated here. Needs to be 1 for all analyses in this group.
reportPCA = 0;                      % Whether PCA stats need to be reported
retinotopyLogic = 'pca';            % Two options here: 'lat' to calculate retinotopy center on response latencies; 'pca' to use early PCA component
reportRetinotopy = 0;               % Whether retinotopy should be reported to the console.
showPCAfigure = 0;                  % Show PCA figure for each brain
showPCAsummaryFigure = 0;           % PCA cumulative figure
showLatDistFigure = 0;              % Show correlations between distance from the center and latency (total figure for all brains)

doEnsembleAnalysis = 0;             % Calcualted adjusted correlations, and identify ensembles from them
doCorrelationFig = 0;               % plot raw correlation matrices - NOT SURE IF UPDATED AT THIS POINT

%%% --- Selectivity group of analyses. The selectivity is always calculated, but we can turn summaries and reporting on and off as we please
showResponseAmplitudes = 0;         % Main simplistic figure for response amplitudes + output of all amplitudes to csv
reportSelectivity = 1;              % Selectivity is always calculated, but this triggers whether it is reported in the console
selectivityName = 'FC';             % Pick which selectivity to report: FC (default), FS, SC, or C (the weighted one, C over [F+S]/2)
showSelTypes = 0;                   % Cell selectivity types figure + output to scv
selToCompare = 'FCtoFS';            % What to compare: FCtoSC (is it a geometry detection?); or FCtoFS (is it a dynamics detection?)
showSelectivityHist = 0;            % Show selectivity histograms (one figure for all brains).

doResponseDynamics = 0;             % Early sells, late cells. Requires doPCA to be on, as it draws from it.
showSpatialEveryBrain = 0;          % A figure for every brain
showSpatialSummary = 0;             % Summary info and summary figure. Note that it doesn't split by age, so run it twice if you need a split

showMugshot = 0;                    % Profile image for every brain

key = 'cfs';
goodSpikeTimeRange = [0.25 2];      % Only this range of times will be kept in spike recordings (to cut out beginning and end artifacts)
msPerFrame = 12;                    % Frames were recorded every 12 ms (about 83 frames/s)



% ------------- Placeholders for global analysis ---------------
averagePcaCurve = [];
goodSpikeTime = [];             % Where good spiking happens. Need to be inialized for consistency.


% ------------- Main cycle -------------
thisPathWithName = mfilename('fullpath');
thisName = mfilename();
localPath = thisPathWithName(1:(length(thisPathWithName)-length(thisName)));    % Path to the folder where this function lies

nBrains = length(folderName);

for(iBrain = 1:nBrains)
    S = [];     % Empty handle for the data structure    
    
    if(age(iBrain)==49); folderBaseIn = 'C:\_Data\___Ca imaging\_caimg s49 mat\';
    else;            folderBaseIn = 'C:\_Data\___Ca imaging\_caimg s46 mat\'; end        
    S = readFiles(folderBaseIn,folderName{iBrain},0);       % Get all raw data in S.f, and all spike reconstructions in S.s. 1st column sweep#, 2nd column time.
    
    nSweeps = length(S);                                % Number of sweeps
    nCells = size(S(1).dataS,2);                        % Number of cells
    time = S(1).timeS;                                  % Good safe time to use later    
    if(isempty(goodSpikeTime))                          % Only build once, to keep consistent across experiments
        goodSpikeTime = find((time>goodSpikeTimeRange(1)) & (time<=goodSpikeTimeRange(2))); % A vector of time points that fall within the good range
    end 
    nTime = length(goodSpikeTime);
    time = time(goodSpikeTime);                         % Only leave middle of the slides, far enough from artifacts.
    timePerTick = mean(time(2:end)-time(1:end-1))*1000; % ms per time tick, used for some images
    if(~randomizePositions)
        xy = S(1).xy;
    else
        xy = S(1).xy(randperm(nCells),:);               % If needed - randomlize everything
        if(iBrain==1); fprintf('NOTE: cell positions were randomized!\n'); end
    end
    % fprintf('%5d cells, %5d stimuli, %5d time points\n',nCells,nSweeps,length(goodSpikeTime));
    if(size(xy,1)>nCells)
        xy = xy(1:nCells,:); 
        S(1).xy = xy; 
        %fprintf('Warning: rows(xy)>nCells in data. Adjusting xy.\n'); 
    end
    name = folderName{iBrain};         % Copy all important info to the output structure    
    
    % ------------- Reformat data for future use (in several different ways) -------------
    fullTrainForOneStimType = [];           % All sweeps spikes for one stim type (pcaStimType), as a 2D matrix, time down fast, sweeps down slow, cells to the right.
    amps = [];                              % Cumulative amplitudes of each response
    ampsPeak = [];                          
    fAvC = [];                              % Averages of fluorescence, across cells
    sAvC = [];                              % Averages of spiking, across cells
    data = [];                              % All sweeps spikes, as a 2D matrix, time down, cells right fast, sweeps right slow
    stimType = mod((1:nSweeps)-1,3)+1;      % Sweep types: 1=C, 2=F, 3=S
    averageShapePerCell = zeros(nTime*3,nCells);  % To keep average traces for every cell, all 3 stim after one another
    for(iSweep=1:nSweeps)
        % figure; plot(S(iSweep).dataF(:,:)); hold on; plot([goodSpikeTime(1) goodSpikeTime(end)],[1 1],'ok'); hold off; error(); % Debugging plot.
        averageShapePerCell((1:nTime)+(stimType(iSweep)-1)*nTime,:) = ...
            averageShapePerCell((1:nTime)+(stimType(iSweep)-1)*nTime,:) + S(iSweep).dataS(goodSpikeTime,:);
        amps = [amps; mean(S(iSweep).dataS(goodSpikeTime,:))];      % Total amplitude of reach cell and each response, nStim by nCells
        fAvC = [fAvC mean(S(iSweep).dataF(goodSpikeTime,:),2)];     % Average fluorescence across entire OT (all cells lumped) for each stimulus
        sAvC = [sAvC mean(S(iSweep).dataS(goodSpikeTime, :),2)];	% Same for spiking
        data = [data S(iSweep).dataS(goodSpikeTime,:)];             % Here concatenation goes right, not down        
    end
    fAvC = bsxfun(@plus,fAvC,-fAvC(1,:));                           % Zero starting points of fluorescence traces
    averageShapePerCell = averageShapePerCell/nSweeps;    
    traceS = repmat(1:3,1,nSweeps/3)';                              % Indicator vector: for each sweep, shows which stimulus type it was.
    
    averageShapePerType = zeros(nTime,3);                           % Average response shapres across entire OT, for this brain
    for(iType=1:3)        
        averageShapePerType(:,iType) = mean(averageShapePerCell((1:nTime)+nTime*(iType-1),:),2);  % Average curves for 3 responses
        [peakY(iType),peakX(iType)] = max(averageShapePerType(:,iType));    % Peak and its position (for plotting and meta-averaging across brains, below)
    end
    
    peakRangePoints = round(peakRange*1000/msPerFrame);               % How far to the left and to the right from the peak to look, for peak estimations
    for(iSweep=1:nSweeps)        
        ampsPeak = [ampsPeak; mean(S(iSweep).dataS(peakX(stimType(iSweep))+(peakRangePoints(1):peakRangePoints(2)),:))]; % Average over a small-ish region around peak
    end    
    
    %figure; % Debugging figure to see whether peak amplitude is 
    %plot(time,averageShapePerType);
    %hold on;
    %plot(time(peakX(2)),0.01,'ks');
    %plot(time(peakX(2))+peakRange,[1 1]*0.01,'k.-');
    %hold off;    
    
    %%% -------------- Show average responses, together with average across all brains --------------
    if(showAverageResponses)
        % Disable averageShapePerType normalization above if you want a picture with correct amplitudes
        flagAverageRspInOnFig = 1;                                              % Do we want them all in one figure, or each in its own figure?
        defaultXShift = 25;                                                     % Default X shift that I just happen to know
        if(flagAverageRspInOnFig)                                               % All in one figure
            if(iBrain==1)
                hF = figure('Color','white','name','avRespPerExp');             % Create figure
                loc.hp(1) = subplot(1,2,1); hold on; xlabel('Time'); ylabel('Responses, s46');    % Axes
                loc.hp(2) = subplot(1,2,2); hold on; xlabel('Time'); ylabel('Responses, s49');                
            else
                hF = findobj('type','figure','name','avRespPerExp');
            end
            
            xShift = peakX(2);              % Peak position foor Flash. We'll try to negate latency variability by shifting the curve by this value
            scalingCoeff = mean(averageShapePerType(:,2));                      % Area for flash resonse
            kAge = (age(iBrain)==49)+1;                                         % Stage code: 1 for stage 46, 2 for stage 49            
            plot(loc.hp(kAge),((1:nTime)-xShift)*timePerTick,averageShapePerType(:,2)/scalingCoeff+00,'-','Color',[1 0.8 0.8]);  % f, 'r'
            plot(loc.hp(kAge),((1:nTime)-xShift)*timePerTick,averageShapePerType(:,3)/scalingCoeff+10,'-','Color',[0.8 1 0.8]);  % s, 'g'
            plot(loc.hp(kAge),((1:nTime)-xShift)*timePerTick,averageShapePerType(:,1)/scalingCoeff+20,'-','Color',[0.8 0.8 1]);  % c, 'b'
            if(iBrain==1)
                averageOfAverages = zeros([size(averageShapePerType) 2]);       % Create empty 3D array with 3d dim representing stages
                avXShift = [0 0];                                               % Average x shifts, for both stages
            end
            for(iStim=1:3)
                shiftedTrace = shift(averageShapePerType(:,iStim)/scalingCoeff,defaultXShift-xShift);           % Shift by that much, padding with zeros (function below)                
                avXShift(kAge) = avXShift(kAge)+xShift;                                                         % Remember this shift, for averaging
               	averageOfAverages(:,iStim,kAge) = averageOfAverages(:,iStim,kAge) + shiftedTrace;
                if(iBrain==nBrains)                                                                             % Last brain; time to average all
                    for(iAge=1:2)
                        thisStage = (iAge==1)*46 + (iAge==2)*49;                                                % Back from stage code to stage (sorry)
                        nBrainsOfThisAge = sum(age(:)==thisStage);                                              % N brains of this age in the set
                        averageOfAverages(:,iStim,iAge) = averageOfAverages(:,iStim,iAge)/nBrainsOfThisAge;     % From sum to average
                        avXShift(iAge) = avXShift(iAge)/nBrainsOfThisAge;
                        yShift = (iStim==1)*20 + (iStim==3)*10;                                                 % How much up to shift the average
                        plot(loc.hp(iAge),((1:nTime)-defaultXShift)*timePerTick,averageOfAverages(:,iStim,iAge)+yShift,'k-');  % Plot                        
                    end
                end
            end
        else                                                                    % Lots of small figures
            figure; hold on; title(name); 
            xShift = 0;                                                         % Plot in original coordinates
            plot(((1:nTime)-xShift)*timePerTick,averageShapePerType(:,1),'b');  % c
            plot(((1:nTime)-xShift)*timePerTick,averageShapePerType(:,2),'r');  % f
            plot(((1:nTime)-xShift)*timePerTick,averageShapePerType(:,3),'g');  % s
        end
        fprintf('%10s\n',name);
        drawnow();
    end
        
    %%% ------------- Normalize response amplitudes within each triad of stimuli ------------
    sampsN = zeros(size(amps));                         % Normalize amplitude responses within each triad
    for(iTriad=1:floor(nSweeps/3))
        ind = (1:3) + (iTriad-1)*3;        
        sampsN(ind,:) = amps(ind,:)/mean(mean(amps(ind,:)));        
        ind = (1:(nCells*3)) + (iTriad-1)*3*nCells;
        k = mean(mean(data(:,ind)));                    % Normalization coefficient
        data(:,ind) = data(:,ind)/k;
        if(pcaStimType>0)
            fullTrainForOneStimType = [fullTrainForOneStimType; ...
                data(:,(1:nCells) + (iTriad-1)*3*nCells + nCells*(pcaStimType-1))];	% Concatenating down, for future PCAing, and in this case Crashes only
        else
            fullTrainForOneStimType = [fullTrainForOneStimType; ...
                data(:,(1:nCells) + (iTriad-1)*3*nCells) ; ...
                data(:,(1:nCells) + (iTriad-1)*3*nCells + nCells) ; ...
                data(:,(1:nCells) + (iTriad-1)*3*nCells + 2*nCells)];           % Concatenating down, for future PCAing, all traces
        end
    end
    
    %%% ------------- Raw data figures -------------
    if(showRawData)
        if(0)
            figure('Color','white');    % --- All sweeps, averaged over cells, per response type
            nsp = ceil(sqrt(nSweeps/3)/1.5);
            msp = ceil(nSweeps/3/nsp);
            for(iTriad=1:floor(nSweeps/3))
                h = axes('units','norm','Position',[floor((iTriad-1)/msp)/nsp mod(iTriad-1,msp)/msp 1/nsp 1/msp]); % Instead of subplot(nsp,msp,iCell), which is extremely slow
                tempc = []; temps = []; tempf = [];
                set(h,'XTickLabel',[],'YTickLabel',[],'nextplot','add');              
                plot( time, sAvC(:,(iTriad-1)*3+1) , 'b-' , 'parent',h);                       
                plot( time, sAvC(:,(iTriad-1)*3+2) , 'r-' , 'parent',h);                       
                plot( time, sAvC(:,(iTriad-1)*3+3) , 'g-' , 'parent',h);                       
            end                 
            drawnow();
        end
        if(0) % Giant one-figure with small colored traces placed in a grid. Not really visible, and takes a while to generate.
            figure('Color','white'); hold on;
            myStyles = 'brg';
            for(iSweep=1:nSweeps)
                for(iCell=1:nCells)
                    plot(time + (iSweep-1)*4, S(iSweep).dataS(goodSpikeTime,iCell)+ (iCell-1)*0.1 , myStyles(mod(iSweep-1,3)+1));
                end 
            end
            title(['All traces; brain ' folderName{iFolder}]);
            drawnow();
        end
        if(1) % Some selected traces, for the paper
            figure('Color','white'); ha1 = axes(); hold on;
            figure('Color','white'); ha2 = axes(); hold on;
            for(q = 1:3)                                % Three sweeps...
                iSweep = [2 3 1]+12;                    % ... but I want them go in this order.
                kx = (time(end)-time(1))*1.1;           % Scaling coeff for x
                for(iCell=1:3)                          % Several traces only
                    ksy = 0.9/max(S(iSweep(q)).dataS(goodSpikeTime,iCell)); % Scaling coeff for y
                    kfy = 0.9/max(S(iSweep(q)).dataF(goodSpikeTime,iCell)-S(iSweep(q)).dataF(goodSpikeTime(1),iCell));
                    plot( ha1, time + (q-1)*kx , S(iSweep(q)).dataS(goodSpikeTime,iCell)*ksy + (iCell-1)*1 , 'b-');
                    plot( ha2, time + (q-1)*kx , (S(iSweep(q)).dataF(goodSpikeTime,iCell)-S(iSweep(q)).dataF(goodSpikeTime(1),iCell))*kfy + (iCell-1)*1 , 'b-');
                end
            end
            hold off;
            drawnow();            
        end
    end    
    
    if(showRawFigure)        
        figure('Color','white');
        triad = 3*((1:floor(nSweeps/3))-1);
        hold on;
        m = mean(sAvC(:,triad+2),2);
        e = std(sAvC(:,triad+2),[],2)/sqrt(nSweeps/3)*tinv(0.025,round(nSweeps/3)-1);
        % Draw a shaded area around the mean curve, in one row (that's why fliplr and everything is a bit weird):
        fill([time' fliplr(time')] , [m'+e' fliplr(m'-e')] , 'r', 'EdgeColor','none', 'FaceAlpha',0.1);
        plot(time, m , 'r-');
        m = mean(sAvC(:,triad+3),2);
        e = std(sAvC(:,triad+3),[],2)/sqrt(nSweeps/3)*tinv(0.025,round(nSweeps/3)-1); 
        fill([time' fliplr(time')] , [m'+e' fliplr(m'-e')] , [0.3 0.7 0], 'EdgeColor','none', 'FaceAlpha',0.1);
        plot(time, m , '-','Color',[0.3 0.7 0]);
        m = mean(sAvC(:,triad+1),2);
        e = std(sAvC(:,triad+1),[],2)/sqrt(nSweeps/3)*tinv(0.025,round(nSweeps/3)-1); 
        fill([time' fliplr(time')] , [m'+e' fliplr(m'-e')] , 'b', 'EdgeColor','none', 'FaceAlpha',0.1);
        plot(time, m , 'b-');        
        hold off;
        xlim([time(1) 1.5]);
        legend({'Flash','','Scramble','','Looming',''});        % Empty '' are used because there are lines and patches on this plat (everything is doubled)
    end
    
    
    %%% ------------- Factor analysis -------------
    if(doPCA)
        m = mean(fullTrainForOneStimType,1);
        fullTrainForOneStimType = fullTrainForOneStimType - m;
        [coeffs,scores,eigs] = pca(fullTrainForOneStimType);        
        scores = scores + repmat(m,size(fullTrainForOneStimType,1),1)*inv(coeffs');
        nComps = 2;    % Either 2 (suggested) or 3
        nEigTo80 = find(cumsum(eigs)/sum(eigs)>=0.8,1,'first'); 
        
        % amps: nStim by nCells
        ampsOfCrashes = amps(stimType==1,:);
        ampsOfCrashes = bsxfun(@plus,ampsOfCrashes,-mean(ampsOfCrashes));
        ampsOfCrashes = bsxfun(@times,ampsOfCrashes,1./std(ampsOfCrashes));
        [~,~,ampEigs] = pca(ampsOfCrashes');                % How different are patterns of total brain response? PCA expects "time" to run right, and presentations - down        
        %figure; subplot(1,2,1); myplot(ampsOfCrashes); subplot(1,2,2); plot(ampEigs/sum(ampEigs)); drawnow(); % Debugging
        nAmpEigTo80 = find(cumsum(ampEigs)/sum(ampEigs)>=0.8,1,'first'); 
        
        if(1) % Rotation - it seems that rotation is always useful
            [~,T] = rotatefactors(scores(:,1:nComps),'Method','promax','Coeff',100);            
            % [~,T] = rotatefactors(coeffs(:,1:nComps),'Method','varimax');
            %T = rotatePositive(scores(:,1:nComps));        % Rotate into quadrant 1. Custom, doesn't work
            coeffs = coeffs(:,1:nComps)*(inv(T)');          % this loooks awkward, but it ensures that s*c'=y still holds (sT(ciT')' = sTiTc'=sc')
            scores = scores(:,1:nComps)*T;                  % Update the scores with promax rotation matrix   
            %T = rotatePositive(coeffs); coeffs = coeffs*T; scores = scores*(inv(T)'); % An attempt to fine-tune coefficients into >0; doesn't work.
            if(sum(scores(:,1))<0); scores(:,1) = -scores(:,1); coeffs(:,1) = -coeffs(:,1); end
            if(sum(scores(:,2))<0); scores(:,2) = -scores(:,2); coeffs(:,2) = -coeffs(:,2); end            
            if(nComps>2); if(sum(scores(:,3))<0); scores(:,3) = -scores(:,3); coeffs(:,3) = -coeffs(:,3); end; end
            if(calculate_latency(scores(:,1))>calculate_latency(scores(:,2)))   % First component should have shorter latency; if it doesn't - swap. Latency is a function defined below.
                scores(:,1:2) = scores(:,1:2)*[0 1; 1 0];
                coeffs(:,1:2) = coeffs(:,1:2)*[0 1; 1 0];
            end
            m1 = median(reshape(scores(:,1),size(scores,1)/nSweeps*3,[]),2);
            m2 = median(reshape(scores(:,2),size(scores,1)/nSweeps*3,[]),2);            
            [~,maxIndm1] = max(m1);             % Positin of max for the first  component
            [~,maxIndm2] = max(m2);             % Positin of max for the second component
            if(maxIndm2<maxIndm1)               % If components are given in a wrong order, rearrange them
                coeffs = coeffs(: , end:-1:1);
                scores = scores(: , end:-1:1);
            end
        end
        
        % -- Report PCA quality results
        if(reportPCA)
            if(iBrain==1)
                fprintf('     Name   Eig1%%   Eig2%%   nEig80  nAmpEigTo80\n');
            end
            fprintf('%10s\t%5.2f\t%5.2f\t%5d\t%5d\n',name,eigs(1)/sum(eigs), eigs(2)/sum(eigs), nEigTo80, nAmpEigTo80);
        end
        
        
        %%% --------- Find the origin point where responses have the shortest latency
        earlyBalance = abs(coeffs(:,1))./(abs(coeffs(:,1))+abs(coeffs(:,2)));      % Response balance in favor of the early component. 
                     % abs() here is to prevent negative weirdos (cells with strange dynamics) from 
                     % creating outliers, as outliers mess with correlations below.
        
        lat = zeros(nCells,1);                          % Place for latencies (note that this code is doubled below in latency analysis, in case this segment isn't run)
        for(iCell=1:nCells)                
            temp = averageShapePerCell((1:nTime)+(pcaStimType-1)*nTime,iCell);
            temp = temp/max(temp);
            lat(iCell) = calculate_latency(temp);                 % Returns a point half-way up to max after 2-piecewise linear fitting                
        end
        lat = lat*msPerFrame;                           % from frames to ms        
                     
        % --- Define function for optimization:
        f = @(x,data) corr(data(:,3),sqrt((data(:,1)-x(1)).^2 + (data(:,2)-x(2)).^2),'Type','Pearson'); % here in data 1-2 will be cell xy, and 3 - its earliness
        opts = optimoptions('fmincon','Display','off');
        switch retinotopyLogic
            case 'pca'
                if(iBrain==1); fprintf('Retinotopy is identified from PCA scores\n'); end;
                originPoint = fmincon(@(a)f(a,[xy earlyBalance(:)]),[60 60],[],[],[],[],[0 0],max(xy,[],1),[],opts);    % Filled args are: starting point, bounds
            case 'lat'
                if(iBrain==1); fprintf('Retinotopy is identified from response latencies\n'); end;                     
                originPoint = fmincon(@(a)f(a,[xy -lat(:)]),[60 60],[],[],[],[],[0 0],max(xy,[],1),[],opts);    % Filled args are: starting point, bounds
                % Optimizing against negative latency as corr with latency with positive, and we are looking for a minimum
        end
        dist = sqrt((xy(:,1)-originPoint(1)).^2 + (xy(:,2)-originPoint(2)).^2);                                 % Distances to the origin point
        %[rDistErl,pDistErl] = corr(dist(:),earlyBalance(:));
        if(reportRetinotopy)
            if(iBrain==1); fprintf('  Name      centX   centY    rDistErl  pDistErl\n'); end;
            fprintf('%8s\t%5.2f\t%5.2f\t%5.2f\t%8s\n',name,originPoint,rDistErl,myst(pDistErl));
        end
        saveResults(auxFolder,folderName{iBrain},'originPoint',originPoint);                                    % Save the output
        %%% position_finder2(xy(:,1),-xy(:,2),earlyBalance); % Good way to visually verify that our "optimal position" actually makes sense        
        
        if(showPCAfigure)
            figure('Color','white','name',name);    % For each brain (not a cumulative figure)
            subplot(2,2,1); hold on;                % --- Response shapes        
            plot(median(reshape(scores(:,1),size(scores,1)/nSweeps*3,[]),2),'b-'); 
            plot(median(reshape(scores(:,2),size(scores,1)/nSweeps*3,[]),2),'r-');
            if(nComps>2); plot(median(reshape(scores(:,3),size(scores,1)/nSweeps*3,[]),2),'g-'); end
            grid on; hold off; %legend({'c1','c2'},'FontSize',6);
            title('1=blue, 2=red');
            subplot(2,2,2); plot(coeffs(:,1),coeffs(:,2),'.'); % --- Loadings
            xlabel('Component 1'); ylabel('Component 2');
            subplot(2,2,3);                         % --- Distribution of cells in space
            if(nComps==2) 
                scatter(xy(:,1),xy(:,2),30,1-[coeffs(:,1) coeffs(:,2) coeffs(:,1)]/max(coeffs(:)),'filled'); set(gca,'Ydir','reverse');
                title('1=green, 2=pink');
            else
                scatter(xy(:,1),xy(:,2),30,1-[coeffs(:,1) coeffs(:,2) coeffs(:,3)]/max(coeffs(:)),'filled'); set(gca,'Ydir','reverse');
                title('1=c, 2=m, 3=y');
            end
            hold on; plot(originPoint(1),originPoint(2),'k+');
            %xlim([min(xy(:,1)) max(xy(:,1))]);        ylim([min(xy(:,2)) max(xy(:,2))]);
            subplot(2,2,4);
            %plot(eigs(1:10)/sum(eigs),'o-'); % Screeplot
            plot(dist,earlyBalance,'.'); % Components in space            
            xlabel('Distance from origin point'); ylabel('early balance');            
            drawnow();
        end
        if(showLatDistFigure)
            if(iBrain==1)
                hF_distlat = figure('Color','White');
                xlabel('Distance'); ylabel('Adjusted Latency, ms');
                ylim([0 2000]);
                hold on;
            else
                figure(hF_distlat);
            end
            plot(dist,lat-mean(lat)+400,'.');
            drawnow();
        end
        
        if(showPCAsummaryFigure)                % Prepare average responses for a summary figure
            [~,maxIndm1] = max(m1);             % Where 1st component maxxes            
            shift = maxIndm1-50;                % How much to shift
            curve1 = median(reshape(scores(:,1),size(scores,1)/nSweeps*3,[]),2);
            curve2 = median(reshape(scores(:,2),size(scores,1)/nSweeps*3,[]),2);
            if(isempty(averagePcaCurve))
                averagePcaCurve = [curve1(:) curve2(:)];
            else
                averagePcaCurve = averagePcaCurve + [curve1(:) curve2(:)];
            end
            if(showRawData)                     % Show two first components for all brains in one figure
                hF = findobj('type','figure','name','allBrains');
                if(isempty(hF)); hF = figure('Color','white','name','allBrains'); hold on; grid on; end;
                m = max(max(curve1(:)-curve1(1)),max(curve2(:)-curve2(1)));
                plot((curve1-curve1(1))/m + (iBrain-1)/2,'b');
                plot((curve2-curve2(1))/m + (iBrain-1)/2,'r');
                drawnow();
            end
        end
    end
    
    
    %%% ----------------- Selectivity calculations ------------
    stimType = mod((1:nSweeps)-1,3)+1;                                      % 1=C, 2=F, 3=S
    itWasaC = stimType==1;
    itWasaF = stimType==2;
    itWasaS = stimType==3;
    if(flagUsePeakAmps)                                                     % Which amplitudes to use: cumulative or peak
        if(iBrain==1)
            fprintf('Warning: using peak amplitudes instead of standard cumulative amplitudes\n');
            fprintf('Looking at the range of %3.1f to %3.1f around the average peak, for each stim type.\n',peakRange);
        end
        amps = ampsPeak;        
    end
    meanRespC = mean(amps(itWasaC,:),1);                                    % "amps" contains total amplitude for every stimulus and every cell (nStim by nCells)
    meanRespF = mean(amps(itWasaF,:),1);                                    % The result of this averaging is a row nCells long
    meanRespS = mean(amps(itWasaS,:),1);
    sdRespC = std(amps(itWasaC,:));
    sdRespF = std(amps(itWasaF,:));
    sdRespS = std(amps(itWasaS,:));
    
    if(showResponseAmplitudes)                                              % Response amplitudes for all cells in the brain.
        if(iBrain==1); giantBagOfAmplitudes = []; end                       % Let's store all mean responses (to 3 stim) of all cells in all brains.
        giantBagOfAmplitudes = [giantBagOfAmplitudes; ...                   % It's saved as a csv later (ctrl-F for the destiny of this variable)
            meanRespF' meanRespS' meanRespC' ones(nCells,1)*iBrain ones(nCells,1)*age(iBrain)]; % f s c ibrain stage
        hF = findobj('type','figure','name','fullAmps');
        if(isempty(hF)); hF = figure('Color','white','name','fullAmps'); hold on; end        
        plot((iBrain-1) + 0.0 + randn(nCells,1)*0.02, meanRespF','r.');
        plot((iBrain-1) + 0.2 + randn(nCells,1)*0.02, meanRespS','.','Color',[0.3 0.7 0]);
        plot((iBrain-1) + 0.4 + randn(nCells,1)*0.02, meanRespC','b.');
        plot((iBrain-1) + [00 0.2 0.4], [mean(meanRespF) mean(meanRespS) mean(meanRespC)],'k-s','MarkerFaceColor','w');
        drawnow;
        if(iBrain==1); fprintf('name       \tnCells      \tnSeeps     \tampF   \tampS   \tampC    \tp(C!=F)    \tp(S!=F)    \tp(C!=S)   \n'); end
        [~,p1] = ttest2(meanRespF,meanRespC);                               % Whether the majority of cells differentiated between C and F
        [~,p2] = ttest2(meanRespF,meanRespS);
        [~,p3] = ttest2(meanRespS,meanRespC);
        fprintf('%10s\t%10d\t%10d\t%10.2f\t%10.2f\t%10.2f\t%10s\t%10s\t%10s\n',name,nCells,nSweeps,...
            mean(meanRespF)*1000,mean(meanRespS)*1000,mean(meanRespC)*1000,myst(p1),myst(p2),myst(p3)); % Multiplied by 1000 to bring it to ~1..10
    end
    
    spiking = mean([meanRespC(:) meanRespF(:) meanRespS(:)],2);             % Most generic average response amplitude (needed for graph_structure analyzer)    
    saveResults(auxFolder,folderName{iBrain},'amp',spiking);                % Save spiking intensity for later analyses.
    selFC = (meanRespC-meanRespF)./sqrt((sdRespC.^2 + sdRespF.^2)/2);       % Cohen d for C-F
    selFS = (meanRespS-meanRespF)./sqrt((sdRespS.^2 + sdRespF.^2)/2);       % same for S-F
    selSC = (meanRespC-meanRespS)./sqrt((sdRespC.^2 + sdRespS.^2)/2);       % same for C-S    
    selC = (meanRespC - (meanRespF+meanRespS)/2)./sqrt((2*sdRespC.^2 + sdRespS.^2 + sdRespF.^2)/4); % Weighted d-like effect size
    
    % selFCp = (mean(ampsPeak(itWasaC==1,:),1)-mean(ampsPeak(itWasaF==1,:),1))./std(ampsPeak(itWasaC==1 | itWasaF==1,:),[],1); % Cohen d for C-F, based on peaks
    
    %%% -- A fancier selectivity (McFadden's Pseudo-R), based on logistic fit
    sel2 = zeros(1,nCells);                             % Needs to be a row for compRow2 calculation below to work.    
    warning('off','all');                               % Glmfit produces annoying warnings when fit is perfect
    for(iCell=1:nCells)
        [~,dev] = glmfit(amps(:,iCell),itWasaC','binomial');   
        devNull = 0;
        nTrialsPseudoR2 = 20;       % How many time to run devNull
        for(q=1:nTrialsPseudoR2)
            [~,temp] = glmfit(amps(:,iCell),itWasaC(randperm(length(itWasaC)))','binomial'); % Repeating this calculation several times actually doesn't help
            devNull = devNull + temp;
        end
        devNull = devNull/nTrialsPseudoR2;
        sel2(iCell) = 1-dev/devNull;        % McFadden's Pseudo-R calculation for logistic regression, and we only take absolute value
    end
    warning('on','all');

    saveResults(auxFolder,folderName{iBrain},'selC',selC);          % Save selectivities for future analysis
    saveResults(auxFolder,folderName{iBrain},'selFC',selFC);
    saveResults(auxFolder,folderName{iBrain},'selSC',selSC);
    saveResults(auxFolder,folderName{iBrain},'selLogistic',sel2);
    
    if(reportSelectivity)        
        switch selectivityName
            case 'FC'; thisSel = selFC(:); fullDelta = mean(meanRespC)/mean(meanRespF)-1;
            case 'FS'; thisSel = selFS(:); fullDelta = mean(meanRespS)/mean(meanRespF)-1;
            case 'SC'; thisSel = selSC(:); fullDelta = mean(meanRespC)/mean(meanRespS)-1;
            case 'C';  thisSel = selC(:);  fullDelta = mean(meanRespC)/((mean(meanRespF)+mean(meanRespS))/2)-1; 
        end
        if(iBrain==1); % Titles in the console
            fprintf('Reporting selectivity: %s\n',selectivityName); 
            fprintf('Experiment       delta   meanSel sdSel  skewness median  90quant    s(>0)  sel2   rSelSel2\n'); 
        end
        fprintf('%12s\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%8.4f\t%5.2f\t%5.2f\n',name,...
            fullDelta, mean(thisSel), std(thisSel), skewness(thisSel,1), median(thisSel), ...
            quantile(thisSel,0.9), sum(thisSel>0)/nCells, mean(sel2), corr(thisSel(:),sel2(:))); 
            %%% Note: Skewness here is not "unbiased" (flag==1), and that's the default for Matlab. It's just a coeff though, not an uncentered moment.
        hF = findobj('type','figure','name','allSel');
        if(isempty(hF)); hF = figure('Color','white','name','allSel'); hold on; xlabel('Brain'); ylabel([selectivityName ' Selectivity']); end
        if(age(iBrain)==49); myColor = 'b.'; else; myColor = 'r.'; end
        plot((iBrain-1) + randn(nCells,1)*0.02,thisSel,myColor);
        plot((iBrain-1),mean(thisSel),'ks','MarkerFaceColor','k');
        drawnow();
    end
    
    if(showSelectivityHist)
        hF = findobj('type','figure','name','selHist');
        if(isempty(hF)); hF = figure('Color','white','name','selHist'); hold on; grid on; end
        figure(hF);
        if(pcaStimType~=3)
            [y,edges] = histcounts(selFC(:)-mean(selFC(:)),-1+(0:30)*0.1);
        else
            [y,edges] = histcounts(selSC(:)-mean(selSC(:)),-1+(0:30)*0.1);
        end
        x = edges(1:end-1)+0.05;                                % Bin centers
        if(iBrain==1); selHist(1:length(x),1) = x; end          % Save bin centers for a summary figure later
        y = y/sum(y);                                           % Normalize
        selHist(1:length(x),iBrain+1) = y;                      % Save data for a summary figure later
        y = y*5+(iBrain-1);                                     % Approximate non-overlap scaling + shift        
        if(age(iBrain)==49); myColor = [0 0 1]; else; myColor = [1 0 0]; end    % Red for 46, blue for 49
        plot(y, x, '.-', 'Color',myColor);
        plot([ones(size(y))*(iBrain-1); y], [x; x], '-','Color',myColor);
        drawnow();        
    end
    
    %%% ------------- Selectivity cell scatters -------------
    if(showSelTypes)        
        hF = findobj('type','figure','name','selTypes');
        if(isempty(hF))
            hF = figure('Color','white','name','selTypes'); 
            ud.a1 = subplot(2,2,1); hold on; xlabel('FS Selectivity');          ylabel('FC selectivity');       set(ud.a1,'FontSize',8);
            ud.a2 = subplot(2,2,2); hold on; xlabel('FS Selectivity adj.');     ylabel('FC selectivity adj');   set(ud.a2,'FontSize',8);
            ud.a3 = subplot(2,2,3); hold on; xlabel('SC Selectivity adj.');     ylabel('FC selectivity adj');   set(ud.a3,'FontSize',8);
            ud.a4 = subplot(2,2,4); hold on; xlabel('SC Selectivity adj.');     ylabel('FC selectivity adj'); 	set(ud.a4,'FontSize',8);
            set(hF,'UserData',ud);            
        else
            ud = get(hF,'UserData');
        end
        if(~exist('giantSelbag','var')); giantSelbag = []; end                  % Let's keep all points to calculate "total r" at the end
        
        [r,p] = corr(selFC(:),selFS(:));
        plot(ud.a1,selFS,selFC,'.');
        plot(ud.a2,selFS-mean(selFS(:)),selFC-mean(selFC(:)),'.');
        plot(ud.a3,selSC,selFC,'.');
        plot(ud.a4,selSC-mean(selSC(:)),selFC-mean(selFC(:)),'.');
        giantSelbag = [giantSelbag; selFC(:) selFS(:) selSC(:) ...          % fc,fs,sc,ibrain,stage
            ones(size(selFC(:)))*iBrain ones(size(selFC(:)))*age(iBrain)];  % Keep all points. See a block governed by 'showSelTypes' somewhere below
        
        if(iBrain==1); fprintf('Name    cor(sel1,sel2)\n'); end
        fprintf('%8s\t%8.2f\t%10s\n',name,r,myst(p));
        
        if(0)                   % Weird figure where responses of different shapes are thrown across the plot based on their selectivity. Doesn't work.
            hold on;
            xscale = 0.01; yscale = 0.06;
            for(iCell=1:nCells)
                ytrace = averageShapePerCell(:,iCell);
                ytrace = ytrace/max(ytrace);
                ytrace = reshape(ytrace,length(ytrace)/3,[]);
                ytrace = ytrace(20:120,:);
                ytrace = resample(ytrace,1,5);
                xtrace = (1:size(ytrace,1));

                %plot(0 + xtrace*xscale , iCell/100+ytrace*yscale,'-');
                plot(selFC(iCell) + xtrace*xscale , selSC(iCell)+ytrace(:,1)*yscale,'-');            
                plot(selFC(iCell) + xtrace*xscale , selSC(iCell)-ytrace(:,3)*yscale,'-');            
            end
        end        
        drawnow();
    end
   
    
     %%% ----------------- Correlation analysis -----------------
    if(doEnsembleAnalysis)    
        corW = zeros(nCells);                           % Future correlation matrix (with average-adjusted correlations)
        for(iStimType=1:3)                              % For every stim type, combine and analyze data
            currentData = [];
            for(iTriad=1:floor(nSweeps/3))      
                currentData = [currentData; ...
                    data(:,(1:nCells) + (iTriad-1)*3*nCells + nCells*(iStimType-1))];	% Concatenating data together            
            end
            tempW = reshuffle_corr(currentData,size(currentData,1)/nTime,0);    % Calculate correlation on small variances
            corW = corW + tempW;                                                % add to the running average
            % if(doCorrelationFig); subplot(3,3,iStimType); myplot(corW); drawnow(); end
        end
        corW = corW/3;                                                          % Finalize average
        
        thresholdType = 'lax';                          % There are several different ways to threshold the corr matrix
        switch(thresholdType)
            case 'lax'
                wThreshold = 0.0;                       % Permissive threshold (only remove negative correlations)
            case 'adaptive'
                mockSet = corW(corW(:)<0);              % Take only negative (presumably false, noise-driven) correlations...
                mockSet = [mockSet(:); -mockSet(:)];    % Make symmetrical, to estimate the width of the "false-positive" curve
                wThreshold = 3*std(mockSet);            % With this threshold, most of fake correlations (assuming that noise is symmetric) will be cut off                
            case 'strict'
                wThreshold = 0.12;                      % At 0.2 everything becomes degenerate and breaks down
        end
        % figure; histogram(corW(:),100); hold on; stem(wThreshold,10,'r'); drawnow(); % Visual check of the threshold value
        
        corW(corW(:)<wThreshold) = 0;               % Remove spurious edges
        if(sum((sum(corW,1)+sum(corW,2)')==0)>0)      % Are there unconnected nodes in corW after thresholding? Answer: no, even with threshold of 0.1.
            warning('W thresholding isolated some nodes from the main graph');
        end
        if(doCorrelationFig)
            figure; subplot(2,3,1); myplot(corW); drawnow();
        end

        [clustInd,modHistory] = spectralClustering(corW,100,10000);            % Clustering analysis (args are maxClusters and sigma parameter)
        
        nClusters = max(clustInd);
        if(doCorrelationFig)             
            [~,ind] = sort(clustInd);
            subplot(2,3,2); myplot(corW(ind,ind)); drawnow(); 
            subplot(2,3,3); hold on;                
            for(i=1:nClusters)
                plot(xy(clustInd==i,1),xy(clustInd==i,2),'o'); 
                %plot(meanRespF(clustInd==i),meanRespC(clustInd==i),'o');
                title(nClusters);
            end
            subplot(2,3,4); plot(modHistory,'.-'); title('Modularity(Nclust)'); drawnow();
        end
        
        %%% Are cells within a cluster closer to each other?
        cell2cellDist = zeros(nCells);          % Distance between any two given cells
        for(ic=1:nCells)
            for(jc=1:nCells)
                cell2cellDist(ic,jc) = sqrt((xy(ic,1)-xy(jc,1))^2 + (xy(ic,2)-xy(jc,2))^2);
            end
        end
        d_within = []; d_between = [];          % Real-world cell-to-cell distances between and within clusters
        nClusters = max(clustInd);
        for(iCluster=2:nClusters)               % Check whether ensembles are closer to each other than between then, in space, and on a graph
            temp = cell2cellDist(clustInd==iCluster,clustInd==iCluster);
            d_within = [d_within; temp(:)];
            for(jCluster=1:(iCluster-1))                
                temp = cell2cellDist(clustInd==iCluster,clustInd==jCluster);
                d_between = [d_between; temp(:)];
                temp = cell2cellDist(clustInd==jCluster,clustInd==iCluster);
                d_between = [d_between; temp(:)];
            end
        end
        clusterCompactness = mean(d_within)/mean(d_between);
        
        if(iBrain==1); fprintf('Name       nClust       maxMod   clustComp\n'); end
        fprintf('%8s\t%5d\t%8.2f\t%8.2f\n',name,nClusters,max(modHistory),clusterCompactness);
        
        saveResults(auxFolder,folderName{iBrain},'corW',corW);
        saveResults(auxFolder,folderName{iBrain},'clustInd',clustInd);
    end
    
    
    
    %%% ------------ Reponse dynamics -----------
    if(doResponseDynamics)
        %%% averageShapePerCell has the following structure: zeros(lgst*3,nCells);
        avResp = mean(averageShapePerCell((1:nTime)+(pcaStimType-1)*nTime,:),2);     % Average response shape across all cells in the tectum, only response to crash
        [~,gloMaxPos] = max(avResp);
        reg1 = (-7:0)+gloMaxPos;                        % 10 ms before the max
        reg2 = (1:8)+gloMaxPos;                         % 10 ms after the max
        reg3 = (9:25)+gloMaxPos;                        % from 10 to 30 ms after the max
        %figure; plot(1:lgst,avResp,'k-'); hold on; plot(gloMaxPos,avResp(gloMaxPos),'o');
        %figure; plot(time,avResp,'k-'); hold on; plot(time(gloMaxPos),avResp(gloMaxPos),'o');
        temp1 = sum(averageShapePerCell(reg1+(pcaStimType-1)*nTime,:),1);     % Amplitude at different time ranges
        temp2 = sum(averageShapePerCell(reg2+(pcaStimType-1)*nTime,:),1);
        temp3 = sum(averageShapePerCell(reg3+(pcaStimType-1)*nTime,:),1);
        aa = [temp1(:) temp2(:) temp3(:)];              % All amplitudes combined
        switch 2
            case 1
                aa = bsxfun(@plus,aa,-min(aa,[],2));                        % Exclude one of the colors
            case 2
                aa = logisticf(bsxfun(@plus,20*aa,-20*mean(aa,2)),0.5,3);   % Increase saturation
        end
        [~,temp] = sort((aa(:,2)-aa(:,1))./(aa(:,2)+aa(:,1)));
        [~,aaRank] = sort(temp);                        % Ranking of color inequality
        %aaRank
        aa = aa/max(aa(:));                             % Norming coeff
        
        % --- Recalculate latency (no guarantee that it was calculated before):
        lat = zeros(nCells,1);
        doLatFig = 0;                               % Troubleshooting set of figures to see whether latency function works        
        for(iCell=1:nCells)
            if(doLatFig); if(mod(iCell-1,20)==0); figure; hold on; end; end
            temp = averageShapePerCell((1:nTime)+(pcaStimType-1)*nTime,iCell);
            temp = temp/max(temp);
            lat(iCell) = calculate_latency(temp);             % Returns a point half-way up to max after 2-piecewise linear fitting
            if(doLatFig); plot(1:length(temp),temp+iCell/2); plot(lat(iCell),temp(lat(iCell))+iCell/2,'ko'); end
        end
        lat = lat*msPerFrame;                       % From frames to ms
        if(doLatFig); drawnow(); end

        [~,temp] = sort(lat);
        [~,latRank] = sort(temp);                       % Ranking by latency
        
        [r1,pval1] = corr(lat(:),selFC(:),'Type','Pearson');               % Correlation: early asymmetry rank to FC selectivity
        [r2,pval2] = corr(dist(:),lat(:),'Type','Pearson');                % Correlation latency vs distance from the origin point
        [r3,pval3] = corr(dist(:),selFC(:),'Type','Pearson');              % Correlation selectivity vs distance from the origin point
        [glmCoeffs,glmDev,glmStats] = glmfit([dist(:) xy(:,1)],selFC(:));
        glmREstimationDist = glmStats.t(2)/sqrt(nCells-2+glmStats.t(2)^2);
        glmREstimationX    = glmStats.t(3)/sqrt(nCells-2+glmStats.t(3)^2);
        if(iBrain==1); fprintf('    Name     rLatSel pLatSel       rLatDist pLatDist     rSelDist pSelDist     x y pGLM estDist estX\n'); end
        fprintf('%10s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%5d\t%5d\t%6s\t%6s\t%6s\n',name,myst(r1),myst(pval1),myst(r2),myst(pval2),myst(r3),myst(pval3),...
            round(originPoint(1)),round(originPoint(2)),myst(glmStats.p(2)),myst(glmREstimationDist),myst(glmREstimationX));
        
        if(showSpatialEveryBrain)   % A figure for every brain
            figure('Color','white');
            subplot(2,2,1);
            top1 = (aa(:,1)>quantile(aa(:,1),1/3));     % Top one-third cells with highest aa1 (early)
            top2 = (aa(:,2)>quantile(aa(:,2),1/3));     % same for mid
            top3 = (aa(:,3)>quantile(aa(:,3),1/3));     % same for late
            hold on;
            plot(((1:nTime)-gloMaxPos) , mean(averageShapePerCell(1:nTime,top1),2),'Color',[0 0.2 1]);
            plot(((1:nTime)-gloMaxPos) , mean(averageShapePerCell(1:nTime,top2),2),'Color',[1 0 0.2]);
            plot(((1:nTime)-gloMaxPos) , mean(averageShapePerCell(1:nTime,top3),2),'Color',[0.5 0.5 0]);
            hold off;
            subplot(2,2,2);
            %scatter(xy(:,1),xy(:,2),30,1-[aa(:,1) aa(:,2) aa(:,3)],'filled'); set(gca,'Ydir','reverse'); title('early=c, mid=m, late=y');
            %scatter(xy(:,1),xy(:,2),30,aaRank(:),'filled'); set(gca,'Ydir','reverse');
            scatter(xy(:,1),xy(:,2),30,latRank,'filled'); set(gca,'Ydir','reverse'); title('early=c, mid=m, late=y');
            subplot(2,3,4);
            plot(latRank(:),selFC(:),'.'); xlabel('Latency rank'); ylabel('Selectivity');
            title(sprintf('r=%s, p=%s',myst(r1),myst(pval1)));
            subplot(2,3,5); 
            plot(dist(:),selFC(:),'.'); xlabel('Distance'); ylabel('Selectivity');
            subplot(2,3,6);            
            %plot(coeffs(:,1)./(coeffs(:,1)+coeffs(:,2)),latRank(:),'.'); % Vs. PCA component balance
            plot(dist,latRank(:),'.'); xlabel('Distance'); ylabel('Latency rank');
            title(sprintf('r=%s, p=%s',myst(r2),myst(pval2)));
            drawnow();
        end
        
        if(showSpatialSummary)  % --------------- Spatial distribution of selectivity
            hF = findobj('type','figure','name','SpatialSummary');            
            if(isempty(hF));                                    % If the figure doesn't exist, create it, and prepare all the axes
                hF = figure('Color','white','name','SpatialSummary'); 
                ud = [];
                ud.ha1 = subplot(1,2,1); hold on; ud.hpMeanLat = plot(1,1,'k.-'); ud.hpBarLat = errorbar(1,1,1,1); 
                xlabel('Latency, ms'); ylabel('Selectivity, d'); ylim([-0.5 Inf]);
                ud.ha2 = subplot(1,2,2); hold on; ud.hpMeanDis = plot(1,1,'k.-'); ud.hpBarDis = errorbar(1,1,1,1); 
                xlabel('Distance, px'); ylabel('Selectivity, d'); ylim([-0.5 Inf]);
                ud.nLat = zeros(1,10); ud.yLat = nan(1,10); ud.vLat = nan(1,10); ud.nDis = zeros(1,10); ud.yDis = nan(1,10); ud.vDis = nan(1,10);
            else
                ud = get(hF,'UserData'); 
            end
            for(q=1:10)
                latBounds = (q-1)*100 + [0 100];                % Latency bounds to bin data, from 0 to 1000
                disBounds = (q-1)*13 + [0 13];                  % Distance, from 0 to 130
                xLat(q) = mean(latBounds);
                yLat(q) = mean(selFC(lat(:)>latBounds(1) & lat(:)<=latBounds(2)));
                xDis(q) = mean(disBounds);
                yDis(q) = mean(selFC(dist(:)>disBounds(1) & dist(:)<=disBounds(2)));
                if(~isnan(yLat(q)))
                    ud.nLat(q) = ud.nLat(q) + 1;                                            % Update the n (number of data point for this average point)                    
                    if(isnan(ud.yLat(q))); ud.yLat(q) = yLat(q); else
                    ud.yLat(q) = (ud.yLat(q)*(ud.nLat(q)-1) + yLat(q))/ud.nLat(q); end      % Update the sum (for mean calculation).                   
                    if(ud.nLat(q)<=2); ud.vLat(q) = 0; else                                 % Note that I'm using a completely messed up way of updating both mean and sd.
                    ud.vLat(q) = (ud.vLat(q)*(ud.nLat(q)-2) + (yLat(q)-ud.yLat(q))^2)/(ud.nLat(q)-1); end % I just got confused that day. It works, but it is 
                end                                                                         % unnecessarily complex, for no reason whatsoever.
                if(~isnan(yDis(q)))
                    ud.nDis(q) = ud.nDis(q) + 1;                                            % Update the n (number of data point for this average point)                    
                    if(isnan(ud.yDis(q))); ud.yDis(q) = yDis(q); else
                    ud.yDis(q) = (ud.yDis(q)*(ud.nDis(q)-1) + yDis(q))/ud.nDis(q); end      % Update the sum (for mean calculation)                    
                    if(ud.nDis(q)<=2); ud.vDis(q) = 0; else
                    ud.vDis(q) = (ud.vDis(q)*(ud.nDis(q)-2) + (yDis(q)-ud.yDis(q))^2)/(ud.nDis(q)-1); end % Update variance
                end   
            end
            ud.xLat = xLat; ud.xDis = xDis;
            set(ud.hpMeanLat,'XData',xLat,'YData',ud.yLat,'Color','k');         % Update mean plot
            ci = sqrt(ud.vLat)./sqrt(ud.nLat).*tinv(0.025,ud.nLat-1);           % Confidence intervals, for error bars
            set(ud.hpBarLat,'XData',xLat,'YData',ud.yLat,'YNegativeDelta',ci,'YPositiveDelta',ci,'Color','k');
            plot(ud.ha1,xLat,yLat,'bo');
            
            set(ud.hpMeanDis,'XData',xDis,'YData',ud.yDis,'Color','k');         % Update mean plot
            ci = sqrt(ud.vDis)./sqrt(ud.nDis).*tinv(0.025,ud.nDis-1);           % Confidence intervals
            set(ud.hpBarDis,'XData',xDis,'YData',ud.yDis,'YNegativeDelta',ci,'YPositiveDelta',ci,'Color','k');
            plot(ud.ha2,xDis,yDis,'bo');
            
            set(hF,'UserData',ud);                              % Upload the data structure
            drawnow;
        end
    end
    
    
    % ------------- Mugshot figure -------------    
    if(showMugshot)
        figure('Color','white');
        set(gcf,'Position',get(gcf,'Position')*[1 0 0 0 ; 0 1 0 -0.3 ; 0 0 1 0 ; 0 0 0 1.2]');
                
        subplot(3,3,1); % --- Fluorescence
        set(gca,'FontSize',8); hold on;
        plot(time,mean(fAvC(:,stimType==1),2),'b-'); % C
        plot(time,mean(fAvC(:,stimType==2),2),'r-'); % F
        plot(time,mean(fAvC(:,stimType==3),2),'g-'); % S
        title(['Date: ' folderName{iBrain}]);        

        subplot(3,3,2); % --- Spiking
        set(gca,'FontSize',8); hold on;
        plot(time,mean(sAvC(:,stimType==1),2),'b-'); 
        plot(time,mean(sAvC(:,stimType==2),2),'r-'); 
        plot(time,mean(sAvC(:,stimType==3),2),'g-'); 
        xlim([time(1) time(end)]);
        title('Spiking');
        legend({'c','f','s'}); legend('boxoff');
        
        %subplot(3,4,3); set(gca,'FontSize',8); % --- History of spiking
        %myplot(amps); colormap('default');
        %title('Amplitude history');
        
        subplot(3,3,3); set(gca,'FontSize',8); % --- History of spiking, adjusted
        [~,ind1] = sort(stimType);                          % To rearrange the amplitudes by stim type
        [~,ind2] = sort(sum(sampsN,1));
        myplot(sampsN(ind1,ind2)); colormap('default');
        title('Adjusted and sorted');

        subplot(3,3,4); set(gca,'FontSize',8); % --- PCA eigenvalues
        [coeffs,scores,eigs] = pca(fullTrainForOneStimType);    % Scores are signals: time down, number right. Coeffs are impacts in each cell (cell down, number right)
        plot(eigs(1:10),'.-');
        title('PCA eigenvalues');
        
        subplot(3,3,5); set(gca,'FontSize',8); % --- FC vs FS
        plot(selFS,selFC,'.'); xlabel('FS selectivity'); ylabel('FC selectivity');
        
        subplot(3,3,6); set(gca,'FontSize',8); % --- FC vs SC
        plot(selSC,selFC,'.'); xlabel('SC selectivity'); ylabel('FC selectivity');

        subplot(3,3,7); set(gca,'FontSize',8); % --- Intensity of spiking in each cell        
        scatter(xy(:,1),xy(:,2),30,mean(fullTrainForOneStimType),'filled');
        xlim([min(xy(:,1)) max(xy(:,1))]);
        ylim([min(xy(:,2)) max(xy(:,2))]);
        title('Intensity');

        %subplot(3,4,6); set(gca,'FontSize',8); % --- THIS IS STUPID AND SHOULD NOT BE USED
        %hold on;
        %plot(mean(reshape(scores(:,2)+scores(:,1)*0.3,lgst*3,[]),2),'b-');
        %plot(mean(reshape(scores(:,3)+scores(:,1)*0.3,lgst*3,[]),2),'r-');
        %title('Guess at response shapes');

        subplot(3,3,8); set(gca,'FontSize',8); % --- Impacts of two components
        rgb = [abs(coeffs(:,3))/max(abs(coeffs(:,3))) zeros(size(coeffs(:,1))) abs(coeffs(:,2))/max(abs(coeffs(:,2)))];
        rgb = 1-repmat(max(rgb,[],2),1,3) + rgb; % Inverses the scale not changing the colors (hard to tell it if it really works)
        scatter(xy(:,1),xy(:,2),30,rgb,'filled');
        xlim([min(xy(:,1)) max(xy(:,1))]);
        ylim([min(xy(:,2)) max(xy(:,2))]);
        title('PCA components 2 and 3');

        subplot(3,3,9); set(gca,'FontSize',8); % --- CS Selectivity
        rgb = [];        
        %S(1).ampS = zeros(nCells,3);                                     % Average full responses to CFS will be stored here. TYPO? ABANDONED?
        for(iCell=1:nCells)
            thisc = (sampsN(mod(traceS,3)==1,iCell));        
            thisf = (sampsN(mod(traceS,3)==2,iCell));
            thiss = (sampsN(mod(traceS,3)==0,iCell));
            %S.ampS(iCell,:) = [mean(thisc) mean(thisf) mean(thiss)];     % S is a vector, so this does not make sense. A TYPO?
            newline = [0 0 0];                                                                      % If nothing is significant, it'll remain black
            if( mean(thisc)>mean(thisf) & ttest(thisc,thisf) ); newline = newline + [0 0 .5]; end   % If C wins, it will be blue.
            if( mean(thisc)>mean(thiss) & ttest(thisc,thiss) ); newline = newline + [0 0 .5]; end   % The more it wins - the bluer it gets.
            if( mean(thiss)>mean(thisc) & ttest(thiss,thisc) ); newline = newline + [0 .5 0]; end   % Winning S makes it green.
            if( mean(thiss)>mean(thisf) & ttest(thiss,thisf) ); newline = newline + [0 .5 0]; end
            if( mean(thisf)>mean(thisc) & ttest(thisf,thisc) ); newline = newline + [.5 0 0]; end   % Winning F - red.
            if( mean(thisf)>mean(thiss) & ttest(thisf,thiss) ); newline = newline + [.5 0 0]; end   
            rgb = [rgb; newline];
        end            
        % traceN = repmat(1:nSweeps,length(time),1);
        % rgb = [mean(spikes(mod(traceN,3)==1,:)); mean(spikes(mod(traceN,3)==2,:)); mean(spikes(mod(traceN,3)==0,:))]';
        % rgb = (rgb-min(rgb(:)))/(max(rgb(:))-min(rgb(:)));
        % rgb = 1-repmat(max(rgb,[],2),1,3) + rgb;            % Inverses the scale not changing the colors (except it doesn't work)
        % rgb = (rgb-repmat(min(rgb),nCells,1))./repmat(max(rgb)-min(rgb),nCells,1); % An attempt to normalize color levels
        % rgb(:,2) = min(rgb(:,1),rgb(:,3));                  % Now R=crash, G=none, B=scrambled
        scatter(xy(:,1),xy(:,2),30,rgb,'filled');
        xlim([min(xy(:,1)) max(xy(:,1))]);
        ylim([min(xy(:,2)) max(xy(:,2))]);
        title('CFS Selectivity');
        drawnow;       
    end
end % file


%%% ---------------------------------------------- SUMMARIES (if any) ----------------------------------

if(showPCAsummaryFigure)
    averagePcaCurve = averagePcaCurve/nBrains;
    figure('Color','white','name','PCA summary'); 
    subplot(1,2,1);
    hold on;
    plot(((1:size(averagePcaCurve,1))-shift)*timePerTick,averagePcaCurve(:,1),'b-'); 
    plot(((1:size(averagePcaCurve,1))-shift)*timePerTick,averagePcaCurve(:,2),'r-');
    legend({'Component 1','Component 2'});
    xlabel('Time, ms'); ylabel('Normalized spiking');
    hold off;
end

if(showSelectivityHist)
    figure;    
    x = selHist(:,1);    y = selHist(:,2:end);
    hold on;
    plot(x,mean(y(:,age==46),2),'r.-');
    plot(x,mean(y(:,age==49),2),'b.-');
    n46 = sum(age==46); 
    n49 = sum(age==49);
    errorbar(x,mean(y(:,age==46),2),std(y(:,age==46),[],2)/sqrt(n46)*tinv(0.025,n46-1),'r.-');
    errorbar(x,mean(y(:,age==49),2),std(y(:,age==49),[],2)/sqrt(n49)*tinv(0.025,n49-1),'b.-');
    hold off;
    xlabel('Selectivity to looming stimulus');
    ylabel('Frequency, %');
    legend({'46','49'});
end

if(showSelTypes)
    %%% ---> expects giantSelbag with 5 columns: FCsel, FSsel, SCsel, brainID, and age
    %figure; subplot(1,2,1); plot(giantSelbag(giantSelbag(:,3)==46,1),giantSelbag(giantSelbag(:,3)==46,2),'.'); title('Stage 46'); 
    %if(selToCompare);     xlabel('FC selectivity'); ylabel('SC selectivity');
    %else; xlabel('FC selectivity'); ylabel('FS selectivity'); end
    %subplot(1,2,2); plot(giantSelbag(giantSelbag(:,3)==49,1),giantSelbag(giantSelbag(:,3)==49,2),'.'); title('Stage 49'); 
    %if(selToCompare);     xlabel('FC selectivity'); ylabel('SC selectivity');
    %else; xlabel('FC selectivity'); ylabel('FS selectivity'); end
    %[r,p] = corr(giantSelbag(:,1),giantSelbag(:,2));
    %fprintf('Total selectivity correlation: %10s\t%s\n',myst(r),myst(p));
    %[r,p] = corr(giantSelbag(giantSelbag(:,3)==46,1),giantSelbag(giantSelbag(:,3)==46,2));
    %fprintf('  s46 selectivity correlation: %10s\t%s\n',myst(r),myst(p));    
    %[r,p] = corr(giantSelbag(giantSelbag(:,3)==49,1),giantSelbag(giantSelbag(:,3)==49,2));
    %fprintf('  s49 selectivity correlation: %10s\t%s\n',myst(r),myst(p));    
    
    fid = fopen([localPath 'sel_allcells_allbrains.csv'],'w');              % csvwrite doesn't support headers, but we need a header
    fprintf(fid,'%s\n','fc,fs,sc,ibrain,stage')                             % so doing it manually
    fclose(fid)
    dlmwrite([localPath 'sel_allcells_allbrains.csv'],giantSelbag,'-append'); % write data to the end
end

if(showResponseAmplitudes)
    %csvwrite([localPath 'avamps_allcells_allbrains.csv'],giantBagOfAmplitudes);
    fid = fopen([localPath 'avamps_allcells_allbrains.csv'],'w');           % csvwrite doesn't support headers, but we need a header
    fprintf(fid,'%s\n','f,s,c,ibrain,stage')                                % so doing it manually
    fclose(fid)
    dlmwrite([localPath 'avamps_allcells_allbrains.csv'],giantBagOfAmplitudes,'-append'); % write data to the end
end

end


function Sbag = readFiles(folderBaseIn,folderName, verbosity)
%%% -------------- Preprocessing -----------------
if(nargin<3); verbosity = 1; end
if(folderName(end)~='\')
    folderName = [folderName '\'];
end
fullFolderName = [folderBaseIn folderName];
if(verbosity); fprintf('\nFolder: %s\n',fullFolderName); end

filesList = dir(fullFolderName);                                        % Find all files in this directory
fileName = [];                                                          % Flush memory
if(length(filesList)>2)
    counter = 0;
    for(q=3:length(filesList))
        if(strcmp(filesList(q).name(end-2:end),'mat'))                  % If matlab data file
            counter = counter + 1;
            fileName{counter} = filesList(q).name;
            % fprintf(' File found: %s\n',fileName{counter});
        end
    end
end
nFiles = length(fileName);

Sbag = [];
for(iFile=1:nFiles)    
    matFileName = [folderBaseIn folderName fileName{iFile}];
    S = [];                                                             % Flush, just in case
    if(verbosity); fprintf('Reading    mat-file: %s\n',matFileName'); end
    load(matFileName,'S');
    if(isfield(S,'key'))                                                % Some S files contain this field, some don't, and it causes problems. 
        S = rmfield(S,'key');                                           % Also I never used it consistenly anyway, so it's not like it is useful for anything.
    end
    try
        Sbag = [Sbag S];
    catch
        S
        Sbag
        error('Concatenation error');
    end
end

end


function T = rotatePositive(scores)
% returns rotation matrix that brings all scores into positive quadrant

k = 1;
[n,m] = size(scores);
T = eye(m);
s = scores;
err0 = sum((s(:)<0)) + k*abs(s(:,1)'*s(:,2));
counter = 0;
figure;

for(q=1:200)
    R = zeros(size(T));
    R(1+floor(length(R(:))*rand(1))) = randn(1)*0.1;    % Only mutate one element in the matrix at a time    
    s = scores*(T+R)/norm(T+R);
    err = sum((s(:)<0)) + k*abs(s(:,1)'*s(:,2));
    if(err<err0)
        T = (T+R)/norm(T+R);
        err0=err;
        counter = counter+1;
        plot(s(:,1),s(:,2),'.');
    end
    xlabel(sprintf('%d %d %5.2f %5.2f',q,counter,err0,err));  drawnow;
end

end


function y = logisticf(x,th,k)
if(nargin<2); th=0.5; end
if(nargin<3); k = 3; end                             % A parameter to set how aggressive the S-shape is. 20 seems to be fine
y = 1./(1+exp((th-x)*k));
% figure; x = 0:0.01:1; plot(x,1./(1+exp((0.5-x)*steepness))); title('sigma-function'); error('just to stop the program');
end


function lat = calculate_latency(trace)
% Calculates onset latency (in frames) for one trace
[~,imax] = max(trace);
imax=imax-5;                            % Move a bit to the left, to make sure we are on the linear segment
f = @(a,x) max(0,a(1)*(x-a(2)));
opt = optimoptions('lsqcurvefit','Display','off');
a = lsqcurvefit(f,[imax/2 2/imax],(1:imax)',trace(1:imax),[1 0],[imax 1000],opt);
lat = round(a(2) + 1/2/a(1));
lat = min(max(1,lat),imax);             % Safety
end


function saveResults(outFolder,brainName,fieldName,value)
% If no file - creates it. If it exists - writes the field.

fullFileName = [outFolder brainName '-auxData.mat'];
try
    D = load(fullFileName);
    D = D.D;                    % I don't like this syntax at all, but load attaches its content as a field, so needs to be de-fielded.
catch
    D = [];
end
D = setfield(D,fieldName,value);
save(fullFileName,'D');
    
end


function data = shift(data,dx)
% Shifts data by dx (integer, positive or negative), padding with zeros.
% Data should be a bunch of vectors with time running down.

if(round(dx)~=dx)
    warning('Can only shift by whole numbers. Rounding now');
    dx = round(dx);
end
[n,m] = size(data);
if(dx>0)
    data = [zeros(dx,m); data(1:end-dx,:)];
elseif(dx<0)
    data = [data(1-dx:end,:); zeros(-dx,m)];
else
    % Do nothing
end

end