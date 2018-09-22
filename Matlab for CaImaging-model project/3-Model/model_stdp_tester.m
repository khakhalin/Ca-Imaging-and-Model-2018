function varargout = model_stdp_tester(type,oneFlag)
% model_stdp_tester
% model_stdp_tester(type,oneFigOnlyFlag)
%
% A tester to study connectivity files created by model_stdp_multisens_1().

% Jun 27 2018: Forked off model_stdp_multisens_1, for dedicated testing. Lots of shared code, unfortunately.
% Jul 03 2018: Now can calculate summaries.
% Jul 18 2018: Daily improvements so far.
% Aug 02 2018: More daily improvements.
% Aug 06 2018: Now functionality for rewiring analysis
% Aug 08 2018: All curves averageing removed (it's now in model_stdp_curve_plotter.m anyway, so no need to duplicate it here)

% Depends on external: 
%   network_rewire (only in nRewires>0 mode, which may be a dead-end), 
%   myst
%   myCentrality - my collection of centrality measures
%   myCyclicity - my attempt to calculate cyclicity
% Currently carries as a copy: myplot


%%% ------------------------------------ Constants ------------------------------------
% rng(3);                                           % Set the random generator seed: keep commented for normal runs

inFolder = 'C:\Users\Arseny\Documents\7_Ca imaging\Model outputs\';  % A folder to look for input files
outFile = ['C:\Users\Arseny\Documents\7_Ca imaging\Model analysis\' datestr(now,'yymmdd') '.csv'];     % A file to put all output values into
if(nargin<1)
    type = 'all';                                   % What type of training files to collect
end
if(nargin<2); oneFlag = 0; end;                     % If not given, assume that we want full analysis, not just one plot
M = [];                                             % In this structure we will keep a table of all data. See functions PUSH, REMEMBER, and the saving block below

nRewires = 1;                                       % If >0, for each normal analysis also performs NREWIRES analyses on randomly rewired data

% The main figure is required as I use it to shuffle handles around, so tolerate it. Can't switch it off.
flagFigDeg = 0;                                     % Whether we want a figure with degree histograms or not
flagFigJunk = 0;                                    % Whether a figure with other junk measures needs to be shown

%%% -------- Read the data
fprintf('Started. Looking for file type: %s\n',type);
fileList = dir(inFolder);                       % List of files
fileList = {fileList(3:end).name};              % Skip '.' and '..' that start the list (at least on Windows)
if(oneFlag); fileList = fileList(1); end;       % If we need only one plot, let there be one plot

%%% -------- Initialize the cumulative figure object
f.hF = figure('color','white');
for(i=1:18)
    f.hp(i) = subplot(3,6,i); set(gca,'FontSize',8); hold on;
end
title(f.hp(1),'Average selectivity');       set(f.hp(1),'YLim',[0 1]);
title(f.hp(2),'Share of sel. cells');       set(f.hp(2),'YLim',[0 1]);
title(f.hp(3),'sum(sel) predict. quality'); set(f.hp(3),'YLim',[0.3 1]);
title(f.hp(4),'Crash Prediction quality');  set(f.hp(4),'YLim',[0.3 1]);
title(f.hp(5),'r Dir-EdgeW');
title(f.hp(6),'Deg dist power');

title(f.hp(7),'%Edge sel2>sel1');                
title(f.hp(8),'2');
title(f.hp(9),'3');
title(f.hp(10),'4');        
title(f.hp(11),'5');

%xlabel(f.hp(13),'Distance');        ylabel(f.hp(13),'Selectivity');
%xlabel(f.hp(14),'Selectivity');     ylabel(f.hp(14),'Influence');
title(f.hp(13),'r Pos-Sel');       xlabel(f.hp(13),'Stage');           set(f.hp(13),'YLim',[-0.5 1]);
%title(f.hp(14),'r Pos-Influ');     xlabel(f.hp(14),'Stage');           set(f.hp(14),'YLim',[-0.5 1]);
title(f.hp(14),'sel90perc');
title(f.hp(15),'bubbliness');      xlabel(f.hp(15),'Stage');

if(flagFigDeg)
    f.hF2 = figure('color','white');
    for(i=1:10)
        f.hp2(i) = subplot(2,5,i); set(gca,'FontSize',8); hold on;
        set(gca,'YScale','log','YLim',[0.01 1],'XLim',[1 8],'XScale','log'); 
    end
end
if(flagFigJunk)
    f.hFj = figure('color','white');
    for(i=1:15)
        f.hpj(i) = subplot(3,5,i); set(gca,'FontSize',8); hold on;
    end
end

f.flagFigDeg = flagFigDeg;
f.flagFigJunk = flagFigJunk;
set(f.hF,'UserData',f);


%%% ------------ Main part
bagResults = [];                                % To keep all OUT values

for(iFile=1:length(fileList))
    try
        temp = load([inFolder fileList{iFile}]);
    catch
        fprintf('Unreadable file: %s\n',fileList{iFile});
        continue;
    end    
    U = temp.U;                                     % Fetch the data    
    if(strcmp(type,U.type) || strcmp(type,'all')) 	% Either correct experiment type or a wildcard
        fprintf('%s  - %s ',fileList{iFile},U.type);
        if(oneFlag)
            M = analyze(M,f.hF,U,5);            % Analyze (main function call), but in this case for one graph only
            close(f.hF);                        % No need to do the summary if there's no summary
            return;          
        else
            for(iStage=1:5)
                tempFileName = fileList{iFile};             % File name
                for(iWire = 1:(1+nRewires))                 % If rewiring is not needed, nRewires==0, and this loop runs only once
                    M = remember(M,'file',tempFileName(1:strfind(tempFileName,'.')-1)); % Cut the extension off
                    M = remember(M,'type',U.type);
                    M = remember(M,'competition',U.competitionMode);
                    M = remember(M,'stage',iStage);
                    if(iWire==1); M = remember(M,'rewire','original'); end % If the graph wasn't rewired, mark it as such
                    if(nRewires>0 && iWire>1);  M = remember(M,'rewire','shuffled'); end % But add it properly
                    M = analyze(M,f.hF,U,iStage,iWire);         % Analyze (main function call)                
                    M = push(M);                                % This row of data is now complete
                    fprintf('.');
                end
            end
        end
        fprintf('\n');
    end
end

%%% -- Time to save the output
hOutFile = fopen(outFile,'w');              % Keep 'w' if we want to overwrite, 'a' to append
for(iCol=1:length(M.title))
    fprintf(hOutFile,'%s',M.title{iCol});
    if(iCol<length(M.title))
        fprintf(hOutFile,',');
    end
end
fprintf(hOutFile,'\r\n');
for(iRow=1:length(M.data))
    fprintf(hOutFile,'%s\r\n',M.data{iRow});
end
fclose(hOutFile);

return;

end



function M = analyze(M,hF,U,iStage,iWire)
% Analyze one of the w matrices (defined by STAGE argument) from this experiment.
% hF - main figure handle
% U - structure from the file, containing w and such
% iStage - which stage to analyze (of 5 stored in the file)
% iWire - if ==1, just analyze it. If >1, rewire for comparison, and mark accordingly.

%%% ------------------ Hard-coded constants ------------------
showFigures = 1;                                % Don't generate figures if called externally (saves time)
nStim = 100;                                    % How many stimuli (of every type) to run
nTickRecalibrate = 2000;                        % Number of ticks (not stim) for network overall excitability (th) recalibration before testing. 1000?
nTick = 15;                                     % Length of one stimulus. For crash, as it is currenly described, 15 is enough

%%% ------------------------------------ Load / reconstruct variables ------------------------------------
% Load some control stuff from the figure:
f = get(hF,'UserData');
% Load model constants from the data variable:
blur = U.blur;
sTarget = U.sTarget;
sAvDecay = U.sAvDecay;
thSens = U.thSens/10;                           % ACHTUNG: That's a rather strong assumption right here! But we need slower change, to make pre-tuning possible.
inactivationMode = U.inactivationMode;

if(iWire==1)
    w = U.stage(iStage).w;
else
    w = network_rewire(U.stage(iStage).w);
end
nCells = size(w,1);
brainDim = sqrt(nCells);
if(brainDim ~= round(brainDim)); error('Non-square brain!'); end;

centerDist = zeros(brainDim,brainDim);                  % Distance from the "virtual tectum" center
for(i=1:brainDim); for(j=1:brainDim); centerDist(i,j) = sqrt((i-brainDim/2-0.5)^2+(j-brainDim/2-0.5)^2); end; end; % Distance from the center of the matrix
edgeDirection = zeros(nCells,nCells);                   % Difference in distance for each potential edge
for(ic=1:nCells); for(jc=1:nCells); edgeDirection(ic,jc) = centerDist(ic)-centerDist(jc); end; end;	
cell2cellDist = zeros(brainDim,brainDim);               % Distance between any two given cells
for(ic=1:nCells)
    ii = floor((ic-1)/brainDim)+1;
    ij = ic-(ii-1)*brainDim;
    for(jc=1:nCells)
        ji = floor((jc-1)/brainDim)+1;
        jj = jc-(ji-1)*brainDim;
        cell2cellDist(ic,jc) = sqrt((ii-ji)^2+(ij-jj)^2);
    end
end

%%% -------------------------------- Recalibration of thresholds before testing -----------------------
% Helps to get rid of short-term effects of most recent learning. Like sleeping overnight. 
% Makes it fair, levels the field for every model, especially after training on different stimuli.
% (Sometimes the most recent stimulus was intense and all neurons are a bit instinsically suppressed, which would bias the test).
% But no hebbian plasticity here; only the thresholds, recalibrated on noise.
sAv = sTarget;                                  % recent average spiking. Assume all spiked ideally at first
th = 1./sTarget/nCells;                         % Starting guess for thresholds
for(q=1:nTickRecalibrate)
    s = floor(rand(nCells,1)+1/nTick);                          % spikes vector (start with something)
    switch inactivationMode
        case 'none'; s = w*s;                                	% Excitatory transmission, no refractory period
        case 'fast'; s = w*s./(1+s);                            % Excitatory transmission + short inactivation period (10 ms)
        case 'slow'; s = w*s./(1+sAv);                          % Excitatory transmission + long inactivation period (~200 ms)
    end
    input = floor(rand(size(s))+2/nTick);                       % Keep the noise at the level of collision stimuli? ACHTUNG: The level of noise is really important!
    s = s + blur*input(:);
    s = logisticf(s,th);                    % Spiking (soft non-linear logistic function, defined below)    
    sAv = sAv*sAvDecay + s*(1-sAvDecay);    % Calculate running average spiking
    th = th - (sTarget-sAv)*thSens;         % Homeostatic plasticity that stays the same
end

%%% ------------------------------------ Testing selectivity ------------------------------------
typeCycle = {'flash','scramble','crash'};       % Not that this sequence of stimuli is different from the sequence used in experiments!
s = sTarget;                                    % spikes vector set to ideal (needs to be reset after training)
sAv = sTarget;                                  % recent average spiking. Assume all spiked ideally at first
nStimTypes = 3;
%amps = zeros(nCells, nTestReps, nStimTypes);    % Total response amplitudes will be stored here
spikeHistory = zeros(nTick,nCells,nStimTypes*nStim);    % Full history of spiking
stimType = zeros(nStimTypes*nStim,1);                   % History of stimuli types
stimCount = 0;
for(si=1:nStimTypes)                                    % For each testing stimulus. 1=F, 2=S, 3=C. Note that ca img experiments go CFS instead, while we go FSC
    inputTrace{si} = zeros(1,nTick);                    % Average input dynamics for this stimulus type (good for troubleshooting)    
    spikeOutput{si} = [];                               % Total spiking (sum over time) of each neuron in each trial
    for(q=1:nStim)                                      % Repeat several times
        stimCount = stimCount+1;
        stimType(stimCount) = si;                       % Remember which stimulus it was
        [mo,~] = generateInput(brainDim,typeCycle{si}); % Create motion object
        s = zeros(size(s));     % Make the network silent
        for(ti=1:nTick)            
            switch inactivationMode
                case 'none'; s = w*s;                               % Excitatory transmission, no refractory period
                case 'fast'; s = w*s./(1+s);                        % Escitatory transmission + short inactivation period (10 ms)
                case 'slow'; s = w*s./(1+sAv);                      % Escitatory transmission + long inactivation period (~200 ms)
            end
            [mo,input] = generateInput(mo);                         % Sensory input
            inputTrace{si}(ti) = inputTrace{si}(ti)+sum(input(:));
            s = s + blur*input(:);
            s = logisticf(s,th);                                        % Spiking (soft non-linear logistic function, defined below)    
            spikeHistory(ti,:,stimCount) = s';                          % Save spiking by each neuron at every moment of time                
        end
    end    
    inputTrace{si}   = inputTrace{si}/nStim;
end

%%% --------- Summaries ------------
%%% spikeHistory : 25 time steps, by nCells (81) cells, by nStimTypes*nStim (75) stimuli
% figure; myplot(squeeze(sum(spikeHistory,1))'); out = 0; return;  % Debugger figure that shows full spiking history

resF = mean(sum(spikeHistory(:,:,stimType==1),1),3);    % Total out of every cell over time, averaged across trials. A vector of numbers, nCell high
resS = mean(sum(spikeHistory(:,:,stimType==2),1),3); 
resC = mean(sum(spikeHistory(:,:,stimType==3),1),3);
spikiness = mean([resF(:) resS(:) resC(:)],2);          % Average demonstrated spikiness; now a column of numbers nCell high

sdF = std(sum(spikeHistory(:,:,stimType==1),1),[],3);   % Standard deviation of cell outputs in each trial, across trials (if resF is MEAN, this one is SD)
sdS = std(sum(spikeHistory(:,:,stimType==2),1),[],3);   % Note that unlike mean(), std() takes dim as argument #3. that's whil this skipped [] argument.
sdC = std(sum(spikeHistory(:,:,stimType==3),1),[],3);

selC = (resC-(resS+resF)/2)./(sqrt((2*sdC.^2+sdS.^2+sdF.^2)/4)); % Weighted d-like effect size
selFC = (resC-resF)./sqrt((sdC.^2+sdF.^2)/2);
selFS = (resS-resF)./sqrt((sdS.^2+sdF.^2)/2);
selSC = (resC-resS)./sqrt((sdC.^2+sdS.^2)/2);


%%% --------- Classification of stimuli by the network ------
responses = squeeze(sum(spikeHistory(:,:,:,1)))';                   % Sum spiking over time, drop dim, flip. Now cells run to the right (cols), stimuli run down (rows)
[loadings,compressedData,~] = pca(responses,'NumComponents',7);     % Reduce dimensionality to make log regression more reliable
actuallyCrash = (stimType==3);                                      % In model files, crash=3, as it goes FSC
firstHalf = repmat([ones(floor(nStim/2),1); zeros(ceil(nStim/2),1)],3,1);                   % First half of responses for each stimulus type
b = mnrfit(compressedData(firstHalf==1,:),actuallyCrash(firstHalf==1)+1,'model','nominal'); % Logistic regression on training half. +1 is needed to turn 0/1 into 1/2 (categories)
prediction = mnrval(b,compressedData(firstHalf~=1,:),'model','nominal');                    % Prediction on the testing half of the dataset
predictionOfCrash = prediction(:,2)>prediction(:,1);
predictionQuality = (sum(predictionOfCrash & actuallyCrash(firstHalf~=1))/sum(actuallyCrash(firstHalf~=1)) + ...
    sum(~predictionOfCrash & ~actuallyCrash(firstHalf~=1))/sum(~actuallyCrash(firstHalf~=1)))/2; % FPR+FNR/2 = expected accuracy on a balanced 50/50 set

%%% Now let's see whether we could detect crashes based on the total tectal output
%%% (Can also be interpreted as "no inhibition", as negative inputs are now not allowed)
bNI = mnrfit(sum(responses(firstHalf==1,selC>0),2),actuallyCrash(firstHalf==1)+1,'model','nominal');    % Logistic regression on training half, and only on selective cells
predictionNI = mnrval(bNI,sum(responses(firstHalf~=1,selC>0),2),'model','nominal');                     % Prediction on the testing half of the dataset
predictionOfCrashNI = predictionNI(:,2)>predictionNI(:,1);
predictionQualityNI = (sum(predictionOfCrashNI & actuallyCrash(firstHalf~=1))/sum(actuallyCrash(firstHalf~=1)) + ...
    sum(~predictionOfCrashNI & ~actuallyCrash(firstHalf~=1))/sum(~actuallyCrash(firstHalf~=1)))/2;

M = remember(M,'fullBrainSel',mean(resC(:))/mean(resF(:))-1);               % Selectivity of the full brain response
M = remember(M,'meanSel',mean(selFC));
M = remember(M,'shareSelCells',mean(selFC>0));
M = remember(M,'sel90perc',quantile(selFC(:),0.9));                         % Threshold for for top 10% selective cells
M = remember(M,'sel90m50',quantile(selFC(:),0.9)-median(selFC(:)));         % Measure of distribution asymmetry (median to 90 gap)
% M = remember(M,'sumOTpredict',predictionQualityNI);
M = remember(M,'bestPredict',predictionQuality);
M = remember(M,'fullBrainSel_SC',mean(resC(:))/mean(resS(:))-1);               % Selectivity of the full brain response
M = remember(M,'meanSel_SC',mean(selSC));
M = remember(M,'shareSelCells_SC',mean(selSC>0));
M = remember(M,'sel90perc_SC',quantile(selSC(:),0.9));                         % Threshold for for top 10% selective cells
M = remember(M,'sel90m50_SC',quantile(selSC(:),0.9)-median(selSC(:)));         % Measure of distribution asymmetry (median to 90 gap)


%%% Note that this logistic regression predicts probability of NOT a crash, as not a crash (0) goes before crash (1).
%%% So negative terms in b mean that respective components are anti-(not-crash), which makes it pro-crash
%%% That's why we need to invert it to get pro-crash prediction power.
%cellInfluence = 1-1./(1+exp(-(b(2:end)'*loadings' + b(1))));   % Transformed to p-like values
cellInfluence = -(b(2:end)'*loadings' + b(1));                  % Raw, untransformed values
cellInfluence = cellInfluence/std(cellInfluence);               % For "crash" the selectivity is so good that it makes the logistic function very sharp, => huge b

rPosSel = corr(centerDist(:),selC(:));                  M = remember(M,'rPosSel',rPosSel);      % Distance from the center COR selectivity
% rPosInf = corr(centerDist(:),cellInfluence(:));         M = remember(M,'rPosInf',rPosInf);      % Distance from the center COR influence
% rSelInf = corr(selC(:),cellInfluence(:));               M = remember(M,'rSelInf',rSelInf);      % Corr selectivition and cell influence (obvious, technical)
rDirWei = corr(edgeDirection(:),w(:));                  M = remember(M,'rDirWei',rDirWei);      % Edge-outwardness (growth of difference from the center) COR its strength
% rDistWei = corr(cell2cellDist(:),w(:));                 M = remember(M,'rDistWei',rDistWei);    % Edge length COR its strength 
mDistWei = sum(sum(w.*cell2cellDist))/sum(w(:))/mean(cell2cellDist(:));                         % Average edge length on our graph, compared to that on a full graph
M = remember(M,'mDistWei',mDistWei);


%%% ------------- Other correlations
% Possible centrality inputs: {'pagerank','revpagerank','netrank','revnetrank','gatherer','reach','revreach','clustering'};
clu = myCentrality(w','clustering');    % Goes from positive for random networks to negative for developed
rSelClu = corr(selC(:),clu);                                M = remember(M,'rSelClu',rSelClu);
rCluSpk = corr(clu(:),spikiness(:));                        M = remember(M,'rCluSpk',rCluSpk);
rSelNet = corr(selC(:),myCentrality(w','netrank'));         M = remember(M,'rSelNet',rSelNet);
rSelRnt = corr(selC(:),myCentrality(w','revnetrank'));      M = remember(M,'rSelRnt',rSelRnt);
rSelGth = corr(selC(:),myCentrality(w','gatherer'));        M = remember(M,'rSelGth',rSelGth);
rSelNIn = corr(selC(:),sum(w,2));                           M = remember(M,'rSelIns',rSelNIn);        % Just a number of inputs to every cell

%selAssort = selC(:)'*w*selC(:)/sum(w(:))-((sum(selC*w)+sum(w*selC'))/(2*sum(w(:))))^2;  % Assortativity: correlation on edges for selectivity - old formula
shi = repmat(selC(:)',nCells,1);                % Matrix of properties for sending neurons
sho = shi';                                     % Matrix of properties for receiving neurons
selAssort = weightedcorr(shi(:),sho(:),w(:));   % Proper weighted correlation
M = remember(M,'selAssort',selAssort);
selComparisonMatrix = bsxfun(@plus,selC(:),-selC(:)');                      % Matrix of sel growth: sel_i - sel_j (to interact with w, where w_ij is an edge from j to i)
shESelGrow = sum((selComparisonMatrix(:)>0).*(w(:)>0.5))/sum(w(:)>0.5);     % What share of strong-ish edges increases Sel?
M = remember(M,'shESelGrow',shESelGrow);
selEGrowth = sum(selComparisonMatrix(:).*w(:))/sum(w(:));                   % Average weighted change in selectivity over an edge
M = remember(M,'selEGrowth',selEGrowth);

%%% --------- Degree distribution analysis -----
temp = sort(w(:));                  % Sort all weights
edgeThreshold = temp(end-nCells);   % find a thershold so that about nCells weights would pass it (similar to how we do it in experiments, to make E=V)
[histDegIn,~] = histcounts(sum(w>edgeThreshold,2),(0:10)-0.5);	% In-degree distributions (of a thresholded graph)
[histDegOu,~] = histcounts(sum(w>edgeThreshold,1),(0:10)-0.5);	% Out-degree distributions
histDegIn = max(histDegIn,1)/nCells;                % Limit above zero to make log() possible
histDegOu = max(histDegOu,1)/nCells;
x = log(1:8);                                       % For fitting, skip deg==0 and deg==1
degFitIn = polyfit(x,log(histDegIn(2:9)),1);        % Fit of in-degrees, from 2 to 8
gammaIn = degFitIn(1);
degFitOu = polyfit(x,log(histDegOu(2:9)),1);        % Fit of out-degrees
gammaOu = degFitOu(1);    
M = remember(M,'gammaIn',gammaIn);
M = remember(M,'gammaOu',gammaOu);
%%% For individual degrees, let's run all analysis on incoming degrees:
deg0  = histDegIn(1)/sum(histDegIn);               M = remember(M,'deg0',deg0);    % Zero inputs
deg12 = sum(histDegIn(2:3))/sum(histDegIn);        M = remember(M,'deg12',deg12);  % 1-2 inputs
deg5p = sum(histDegIn(6:end))/sum(histDegIn);      M = remember(M,'deg5p',deg5p);  % more than 5 inputs

if(f.flagFigDeg)    
    plot(f.hp2(0+iStage),0:9,histDegIn,'kx');    
    plot(f.hp2(0+iStage),1:8,exp(degFitIn(2))*(1:8).^(gammaIn),'-','Color',[1 0.9 0.9]);    
    plot(f.hp2(5+iStage),0:9,histDegOu,'kx');
    plot(f.hp2(5+iStage),1:8,exp(degFitOu(2))*(1:8).^(gammaOu),'-','Color',[1 0.9 0.9]);
    %plot(f.hp2(5+iStage),2:8,c2(2:8),'-','Color',[0.9 0.9 1]);
end

%myFitType = fittype('exp(a*log(x))','independent','x');                % Sense-check. Fitting in log-space yields slightly different results than fitting in raw p-k space.
%myFitOptions = fitoptions('Method','NonlinearLeastSquares','Lower',[-100],'Upper', [10],'Startpoint',[-1]); % But that's OK in this situation, I think.
%[c2,gof2] = fit((2:8)',max(histDegOu(3:9)',1)/nCells,myFitType,myFitOptions);
%c2
%fprintf('Stage %d, gammaIn %f, gammaOu %f\n',iStage,gammaIn,gammaOu); % Troubleshooting


%%% --------- Dynamic richness of responses -----
% The question here: are all responses of the same type alike?
% Cannot do on "S" though, as in the current implementation they are truly trial-by-trial randomized (unlike in the experiment)
% Key variables: spikeHistory(time#,neuron#,stim#), stimType (a vector of 123s)
respF = zeros(nTick*sum(stimType==1),nCells);
respC = zeros(nTick*sum(stimType==3),nCells);
counter = 0;
for(iStim=find(stimType'==1))    % Go through all F stimuli
    respF(counter*nTick+(1:nTick),:) = squeeze(spikeHistory(:,:,iStim));
    counter = counter+1;
end
counter = 0;
for(iStim=find(stimType'==3))    % Go through all C stimuli
    respC(counter*nTick+(1:nTick),:) = squeeze(spikeHistory(:,:,iStim));
    counter = counter+1;
end
[~,~,eigenV] = pca(respF); nRichTo80F = find(cumsum(eigenV)/sum(eigenV)>=0.8,1,'first'); M = remember(M,'nRichTo80F',nRichTo80F);
[~,~,eigenV] = pca(respC); nRichTo80C = find(cumsum(eigenV)/sum(eigenV)>=0.8,1,'first'); M = remember(M,'nRichTo80C',nRichTo80C);

%%% --------- Ensembles -----
% Key variables: spikeHistory(time#,neuron#,stim#), stimType (a vector of 123s)
% We should start with preparing unbiased (centered) and normalized contribution values.
ampF = squeeze(sum(spikeHistory(:,:,stimType==1),1));   % Total output of each cell over time, in each trial: Neurons down, stimuli right
ampS = squeeze(sum(spikeHistory(:,:,stimType==2),1));
ampC = squeeze(sum(spikeHistory(:,:,stimType==3),1));
ampF = bsxfun(@plus,ampF,-mean(ampF,1)); ampF = bsxfun(@times,ampF,1./std(ampF,[],1));  % Unbiased an normalized
ampS = bsxfun(@plus,ampS,-mean(ampS,1)); ampS = bsxfun(@times,ampS,1./std(ampS,[],1));
ampC = bsxfun(@plus,ampC,-mean(ampC,1)); ampC = bsxfun(@times,ampC,1./std(ampC,[],1));
%[loadings,compressedData,~] = pca(ampC,'NumComponents',10);  
[~,compressedData,eigenV] = pca([ampF ampS ampC]);  
compressedData = compressedData(:,1:10);                        % For clustering, save only 10 first compoments
nPCAto80 = find(cumsum(eigenV)/sum(eigenV)>=0.8,1,'first');     % How many components one would need to explain 80% of variance
M = remember(M,'nPCAto80',nPCAto80);

if(1)
    [clustInd,varExplained] = kmeans_opt(compressedData);                   % Attempt of optimal clustering        
else
    [~,~,fullVar]      =kmeans(compressedData,1,'emptyaction','drop');      % Alternatively - go with an arbitrary "reasonable" number of clusters (not recommended)
    [clustInd,~,resVar]=kmeans(compressedData,8,'emptyaction','drop');      % Default distance is "squared Euclidean", so no need to square
    varExplained = 1-sum(resVar)/fullVar;
end
% figure; scatter(compressedData(:,1),compressedData(:,2),10,clustInd,'filled'); % Visual debugging - cluster-gram.
% [~,ind] = sort(clustInd);
% figure; myplot(ampC(ind,:));       % Visual debugging. Will we see the bands?
% myCor = corr([ampF ampS ampC]');
% figure; myplot(myCor(ind,ind));    % Visual debugging. Will we see the quilt?
w_within = []; w_between = [];          % Weight of edges connecting within and between clusters
d_within = []; d_between = [];          % Real-world cell-to-cell distances between and within clusters
nClusters = max(clustInd);
for(iCluster=2:nClusters)               % Check whether ensembles are closer to each other than between then, in space, and on a graph
    temp = w(clustInd==iCluster,clustInd==iCluster);
    w_within = [w_within; temp(:)];
    temp = cell2cellDist(clustInd==iCluster,clustInd==iCluster);
    d_within = [d_within; temp(:)];
    for(jCluster=1:(iCluster-1))
        temp = w(clustInd==iCluster,clustInd==jCluster);
        w_between = [w_between; temp(:)];
        temp = w(clustInd==jCluster,clustInd==iCluster);
        w_between = [w_between; temp(:)];
        temp = cell2cellDist(clustInd==iCluster,clustInd==jCluster);
        d_between = [d_between; temp(:)];
        temp = cell2cellDist(clustInd==jCluster,clustInd==iCluster);
        d_between = [d_between; temp(:)];
    end
end    
[~,pval] = ttest2(w_within,w_between);
% fprintf('N clusters: %2d; Within: %5.3f; between: %5.3f (p=%s); var explained: %4.2f\n',...
%     nClusters,mean(bagWithin),mean(bagBetween),myst(pval),varExplained); % Console reporting
withinClusterPreference = mean(w_within)/mean(w_between);
clusterCompactness = mean(d_within)/mean(d_between);
M = remember(M,'nClusters',nClusters);
M = remember(M,'clustVarExplained',varExplained);
M = remember(M,'clusterPreference',withinClusterPreference);
M = remember(M,'clusterCompactness',clusterCompactness);
M = remember(M,'clustPrefPval',pval);

%%% --------- Other misc measures ------
% Key: 1 out-in, 2 in-out, 3 out-out, 4 in-in
% asOI = assortativity_wei(U.stage(iStage).w',1); M = remember(M,'asOI',asOI);
% asIO = assortativity_wei(U.stage(iStage).w',2); M = remember(M,'asIO',asIO);
% asOO = assortativity_wei(U.stage(iStage).w',3); M = remember(M,'asOO',asOO);
% asII = assortativity_wei(U.stage(iStage).w',4); M = remember(M,'asII',asII);

%asOI = myCentrality(U.stage(iStage).w','assoi');    M = remember(M,'asOI',asOI);
%asIO = myCentrality(U.stage(iStage).w','assio');    M = remember(M,'asIO',asIO);
%asOO = myCentrality(U.stage(iStage).w','assoo');    M = remember(M,'asOO',asOO);
%asII = myCentrality(U.stage(iStage).w','assii');    M = remember(M,'asII',asII);

eff = efficiency_wei(U.stage(iStage).w');           M = remember(M,'eff',eff);
[~,modul] = modularity_dir(w');                     M = remember(M,'modul',modul);
temp = clustering_coef_wd(w'); temp = mean(temp(~isinf(temp)));  M = remember(M,'clust',temp);
[flow, revFlow] = network_flow(w');                 M = remember(M,'flow',flow);
                                                    M = remember(M,'revFlow',revFlow);
cycl = myCyclicity(U.stage(iStage).w');             M = remember(M,'cycl',cycl);

rSelSpk = corr(selC(:),spikiness(:));               M = remember(M,'rSelSpk',rSelSpk);
rSelfcSelfs = corr(selFC(:),selFS(:));              M = remember(M,'rSelfcSelfs',rSelfcSelfs);
rSelfcSelsc = corr(selFC(:),selSC(:));              M = remember(M,'rSelfcSelsc',rSelfcSelsc);

%%% --------- Visualize outputs -----
markerStyleList = {'bo','ro','kx','gs','mx'};   % A piece of code here is to cycle through marker styles when vizualizing a mix of different renders
switch U.competitionMode % in, out, global, decay, none, softout, satin
    case 'in';      markerStyle = 'b.';
    case 'out';     markerStyle = 'r.';
    case 'global';  markerStyle = 'kx';
    case 'softin';  markerStyle = 'bo';
    case 'softout'; markerStyle = 'ro';
    case 'satin';   markerStyle = 'ms';
    case 'satout';  markerStyle = 'gs';
    case 'slide';   markerStyle = 'bx';
    otherwise;      markerStyle = 'k.';
end
if(iWire>1)     % Reshuffled analysis - override marker style
    markerStyle = 'c.';
end

plot(f.hp(1),iStage,mean(selC),markerStyle);            
plot(f.hp(2),iStage,mean(selC>0),markerStyle);                                                             
plot(f.hp(3),iStage,predictionQualityNI,markerStyle);   
plot(f.hp(4),iStage,predictionQuality,markerStyle);     
plot(f.hp(5),iStage,rDirWei,markerStyle);               % already saved above
plot(f.hp(6),iStage,(gammaIn+gammaOu)/2,markerStyle);   % saved above, in and out separately
plot(f.hp(7),iStage,shESelGrow,markerStyle);            % saved above
plot(f.hp(8),iStage,selEGrowth,markerStyle);            title(f.hp(8),'selGrowth');
plot(f.hp(9),iStage,nRichTo80F,markerStyle);            title(f.hp(9),'Dyn Rich F');
plot(f.hp(10),iStage,nRichTo80C,markerStyle);           title(f.hp(10),'Dyn Rich C');
plot(f.hp(11),iStage-0.2,deg0,'b.');    plot(f.hp(11),iStage,deg12,'k.');   plot(f.hp(11),iStage+0.2,deg5p,'r.');	title(f.hp(11),'Degrees');
plot(f.hp(13),iStage,rPosSel,markerStyle);                  % All thesea are already saved above
plot(f.hp(14),iStage,quantile(selC(:),0.9),markerStyle);    % 90s percentile of selectivity
plot(f.hp(15),iStage,withinClusterPreference,markerStyle);
plot(f.hp(16),iStage,rSelClu,markerStyle); title(f.hp(16),'r Sel Clust');       
plot(f.hp(17),iStage,rSelNet,markerStyle); title(f.hp(17),'r Sel NetRank');     
plot(f.hp(18),iStage,rSelRnt,markerStyle); title(f.hp(18),'r Sel RevNetRank');  

%%% Excluded plots:
%plot(f.hp(14),iStage,rPosInf,markerStyle);
%plot(f.hp(15),iStage,rSelInf,markerStyle);

if(f.flagFigJunk)
    plot(f.hpj(1),iStage,asOI,markerStyle);         title(f.hpj(1),'out-in');   % These vars are all M-saved above
    plot(f.hpj(2),iStage,asIO,markerStyle);         title(f.hpj(2),'in-out');
    plot(f.hpj(3),iStage,asOO,markerStyle);         title(f.hpj(3),'out-out');
    plot(f.hpj(4),iStage,asII,markerStyle);         title(f.hpj(4),'in-in');    
    plot(f.hpj(5),iStage,nPCAto80+randn(1)*0.1,markerStyle);    title(f.hpj(5),'nPCAto80');    % Some jitter is necessary here    
    plot(f.hpj(6),iStage,cycl,markerStyle);         title(f.hpj(6),'cyclic'); set(f.hpj(6),'YScale','log');
    plot(f.hpj(7),iStage,rSelSpk,markerStyle);      title(f.hpj(7),'rSelSpk');
    plot(f.hpj(8),iStage,rCluSpk,markerStyle);      title(f.hpj(8),'rCluSpk');
    plot(f.hpj(9),iStage,selAssort,markerStyle);    title(f.hpj(9),'selAssort');  
    plot(f.hpj(10),iStage,rSelGth,markerStyle);     title(f.hpj(10),'r Sel Gather');    
    plot(f.hpj(11),iStage,rDistWei,markerStyle);    title(f.hpj(11),'rDistWei');
    plot(f.hpj(12),iStage,mDistWei,markerStyle);    title(f.hpj(12),'dConn/dAll');
    %plot(f.hpj(13),iStage,rSelNIn,markerStyle);     title(f.hpj(13),'rSelIns');
    plot(f.hpj(13),iStage,eff,markerStyle);         title(f.hpj(13),'Efficiency');
    plot(f.hpj(14),iStage,clusterCompactness,markerStyle);     title(f.hpj(14),'clust compact');
    plot(f.hpj(15),iStage,rSelfcSelsc,markerStyle); title(f.hpj(15),'rSelfcSelsc');
end

if(iStage==5)    % Well-developed network
    %plot(f.hp(13),centerDist(:),selC,'b.');    
    %plot(f.hp(14),selC,cellInfluence(:),'b.');    
    
    %scatter(f.hp(15),sum(w,2),sum(w,1)',5,selC,'filled');
    if(0) % Various one-time figures
        %figure; scatter(centerDist(:),sum(w,1),10,selC,'filled'); % Position, sum of outputs, and selectivity
        %figure; scatter(centerDist(:),selC,15,sum(w,1),'filled'); % Position, sum of outputs, and selectivity
        
        selectivity_graph(w,selC); % <<------------------------- Call to another huge chunk of analysis
    end
    
    if(0) % Troubleshooting segment
        [loadings selC(:) cellInfluence(:)]
        b
        figure;
        for(i=1:7)
            subplot(2,4,i); myplot(reshape(loadings(:,i),brainDim,brainDim));
            title(sprintf('b = %4.2f',b(i+1)));
        end
        subplot(2,4,8); myplot(reshape(cellInfluence,brainDim,brainDim));
    end
end
%plot(f.hp(8),compressedData(firstHalf~=1,:)*b(2:end),prediction(:,1),'b.');
drawnow(); 


%%% --------- Send same outputs out ---------
return;

%%% -------------------------------------------------------------------- end of it for now, the rest is abandoned --------------------






if(showFigures) % --- How selectivity interacts with cell place within the network
    figure('Color','white');
    lineStyle = {'k-','r-','b-'};
    %h6 = subplot(3,3,6); hold on;
    for(si=1:length(spikeHistory))
        subplot(3,3,si);   myplot(reshape(mean(spikeOutput{si},2),dim,dim)); title(typeCycle{si}); colorbar(); %caxis([0 2]);
        subplot(3,3,3+si); myplot(spikeHistory{si}); title(typeCycle{si}); %caxis([0 2]);    
        %plot(h6,inputTrace{si},lineStyle{si});    
    end
    subplot(3,3,6); plot(resC-resF,selFC,'.'); xlabel('C minus F'); ylabel('Selectivity');
    subplot(3,3,3); hold on; plot([0 2],[0 2],'g-'); plot(resF,resC,'.'); hold off; xlabel('to flash'); ylabel('to crash'); title('Cell responses');
    subplot(3,3,7); plot(sum(w), selFC,'.'); xlabel('total syn drive'); ylabel('Selectivity');
    subplot(3,3,8); plot(sTarget,selFC,'.'); xlabel('spike target');    ylabel('Selectivity');
    subplot(3,3,9); plot(th,     selFC,'.'); xlabel('spike threshold'); ylabel('Selectivity');
    drawnow();
end

if(0)   % Long console output
    fprintf('Share of neurons selective to crash: %f\n',sum(resC>resF)/nCells);    
    fprintf('Median crash/flash: %f\n',median(resC./resF));
    fprintf('Average crash minus flash: %f\n',mean(resC-resF));    
end
if(1)   % Short console output
    if(showFigures) % The reason this is linked to show figures is that it is ==0 during serial runs of the model
        fprintf('nStim, nCollisions, share Collisions, Full spiking F, full spiking C, share of C>F, share of C>120F, average selFC, var selFC, dist|sel\n');
    end    
    distanceToCenter = sqrt((meshx(:)-dim/2).^2 + (meshy(:)-dim/2).^2);
    [rho_dsel,pval_dsel] = corr(distanceToCenter,selFC);
    fprintf('%5d\t%5d\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%5.2f', nStim , nCollisions , sum(resF) , sum(resC),...
        sum(resC>resF)/nCells , mean(selFC) , var(selFC), median(selFC), quantile(selFC,0.9), rho_dsel);
    %for(q=1:length(newDataLine))
    %    fprintf('\t%7f',newDataLine(q));
    %end
    fprintf('\n');
end

if(nargout>0); varargout{1} = w; end;
if(nargout>1); varargout{2} = selFC; end;
if(nargout>2); varargout{3} = th; end;

if(showFigures) % --- Show final graph in both original and optimized coordinates, 
                % as well as whether selectivity correlates with any centrality measure
    selectivity_graph(w,selFC,meshx,meshy);   % Analyze the graph and show some more figures (see the procedure itself)
    
    figure; % Selectivity as a function of distance    
    [rho,pval] = corr(distanceToCenter,selFC);
    plot(distanceToCenter,selFC,'.'); xlabel('Distance from the center'); ylabel('FC Selectivity');
    title(sprintf('r=%s; p=%s',myst(rho),myst(pval)));
    
    figure; hold on; % Alternative figure for full responses
    plot([1 2 3],[resF(:) resS(:) resC(:)],'-','Color',[1 1 1]*0.9);
    plot(1,resF,'r.');  plot(1,mean(resF),'ks');
    plot(2,resS,'g.');  plot(2,mean(resS),'ks');
    plot(3,resC,'b.');  plot(3,mean(resC),'ks');  
    hold off; xlabel('Stimulus type'); ylabel('Cell response'); xlim([0 4]);
end

end



function y = logisticf(x,th)
steepness = 20;                             % A parameter to set how aggressive the S-shape is. 20 seems to be fine
y = 1./(1+exp((th-x)*steepness));
% figure; x = 0:0.01:1; plot(x,1./(1+exp((0.5-x)*steepness))); title('sigma-function'); error('just to stop the program');
end


function [mo,s] = generateInput(varargin)
% [mo,~] = generateInput(dim,type) - generates a motion settings object
% [mo,s] = generateInput(mo) - generates an inpus sequence S, and updates the motion object

% rgcDecay = 0;  % Real RGCs burst for about 200 ms (until response decreases 10 times) k = exp(log(0.1)/20) = 0.9
                 % In fact decay totally prevents selectivity from happening. The topology is similar, but no selectivity.
arguments = varargin(1:end);
if(length(arguments)>1) % --- Initialize: create a motion object
    dim = varargin{1};
    mo.type = varargin{2};       
    mo.r = dim/4;                                       % Circle radius
    mo.period = 10;                                     % Collision length. 10 ticks is 100 ms    
    mo.shuffle = 0;                                     % Default value: no shuffle
    switch mo.type
        case 'flash'
            mo.dstart = 0.4;                                    % Starting distance to the object
            mo.dfinal = -3;                                     % Final distance (doesn't make much sense in this case). Adjusted to make the transition last 2 frames
            mo.startAngle = 0;                                  % Starting angle
            mo.moveAngle = 0;                                   % Moving direction (opposite to the starting angle + noise)
            mo.start = [dim/2, dim/2]+0.5;                      % Starting position
            mo.v = 0;                                           % Lateral speed
            mo.shuffle = 1;                                     % This one needs to be shuffled, to have some variability in responses
        case 'crash'
            mo.dstart = 1;                                      % Starting distance to the object.
            mo.dfinal = 0.0001;                                 % Finishing distance to the object (never 0). Doesn't really matter, as it's overridden below
            mo.startAngle = 0;                                  % Starting angle
            mo.moveAngle = 0;                                   % Moving direction (opposite to the starting angle + noise)
            mo.start = [dim/2 dim/2]+0.5+randn(1,2);            % Starting position
            mo.v = 0;                                           % Lateral speed
        case 'scramble'
            mo.dstart = 1;                                      % Starting distance to the object.
            mo.dfinal = 0.0001;                                 % Finishing distance to the object (never 0). Doesn't really matter, as it's overridden below
            mo.startAngle = 0;                                  % Starting angle
            mo.moveAngle = 0;                                   % Moving direction (opposite to the starting angle + noise)
            mo.start = [dim/2 dim/2]+0.5+randn(1,2);            % Starting position
            mo.v = 0;                                           % Lateral speed
            mo.shuffle = 1;                                     % Shuffle this one
        case 'trans' % Right now unbalanced in terms of total retinal input
            mo.dstart = 0.5;                                    % Starting distance to the object
            mo.dfinal = 0.5;                                    % Finishing distance to the object (never 0)            
            mo.startAngle = pi;                                 % Starting angle
            mo.moveAngle = 0;                                   % Moving direction (opposite to the starting angle + noise)
            mo.start = [0,dim/2]+0.5;                           % Starting position   
            mo.v = dim/mo.period;                               % Speed. It will cross the field in this time
        case {'vis','shuffle','full'} % Training transitions
            mo.dstart = 1;                                      % Starting distance to the object
            mo.dfinal = rand(1)*1+0.01;                         % Finishing distance to the object (never 0)    
            mo.startAngle = rand(1)*2*pi;                       % Starting angle
            mo.moveAngle = mo.startAngle+pi+rand(1)*pi/3-pi/6;  % Moving direction (opposite to the starting angle + noise)
            offset = rand(1)*(dim/2+mo.r-1);                    % Starting offset from the center (from 0 to not visible)
            mo.start = [dim/2 + cos(mo.startAngle)*offset, dim/2 + sin(mo.startAngle)*offset]+0.5;  % Starting position
            mo.v = dim/mo.period;                               % Speed. It will cross the field in this time
            if(strcmp(mo.type,'shuffle'))
                mo.shuffle = 1;                                 % Set Shuffle to 1 if necessary
            end
        case 'rand' % Just fill with random noise
            mo.dstart = 1; mo.dfinal = 1; mo.startAngle = 0; 
            mo.moveAngle = 0; mo.start = [0 0]; mo.v = 0; 
    end
    
    mo.d = mo.dstart;                                   % Current distance
    mo.dstep = (mo.dfinal-mo.dstart)/mo.period;         % Distance change speed
    mo.pos = mo.start;
    mo.dir = [cos(mo.moveAngle) sin(mo.moveAngle)];     % Moving vector (for now - opposite of original location)        
    
    [mo.i,mo.j] = meshgrid(1:dim,1:dim);                % A technical variable, used during actual drawing
    mo.retina = zeros(dim);                             % What the retina sees
    mo.rgc = zeros(dim);                                % Output of RGCs (a change in retina + bursting)    
    mo.collided = 0;                                    % Ever collided (for this object)
    mo.justCollided = 0;                                % Collided this frame
    mo.emptyFrames = 0;                                 % Empty frames in a row (to save time if the object left)
    mo.noise = 1/dim^2;                                 % If non-zero, adds XOR noise to the signal. 1/dim^2 means on average one noisy pixel per frame
    s = zeros(dim);                                     % Even though it's not used at init, populate with something to prevent from crashing
    
else  % --- Generate new frame using motion object "mo"
    mo = varargin{1};
    nCells = numel(mo.retina);
    dim = sqrt(nCells);    
    
    switch mo.type
        case 'rand'
            noiseLevel = 1/dim;      % Noise level.
            s = floor(rand(nCells,1)+noiseLevel);
        otherwise         
            if(~mo.collided)                            % If not collided yet
                oldretina = mo.retina;                                      % Save old retina
                mo.retina = zeros(dim);                                     % Start with a white board
                switch mo.type
                    case {'crash','scramble'}
                        currentR = sqrt(2)*dim/2*(1-mo.d);                  % Linear transition from 0 to full screen
                    otherwise
                        currentR = mo.r/mo.d;                               % Realistic radius (perceived) right now
                end
                mo.retina((mo.i-mo.pos(1)).^2 + (mo.j-mo.pos(2)).^2 <= currentR^2) = 1; % Calculate what is seen (draw a circle)
                mo.d = mo.d+mo.dstep;                                       % Change distance to the object
                mo.pos = mo.pos + mo.dir*mo.v;                              % Move the circle to a new place
                mo.rgc = xor(mo.retina,oldretina);                          % Only report change (on/off cells)
                % mo.rgc = mo.rgc*rgcDecay + xor(mo.retina,oldretina);      % Older version with decay; doesn't work (see comments)                
                s = mo.rgc;
                if(strcmp(mo.type,'full') || strcmp(mo.type,'vis') || strcmp(mo.type,'shuffle')) % Training stimuli, not testing
                    sumRet = sum(mo.retina(:));
                    if(sumRet>0.75*nCells & mo.collided==0)                 % Collision happened (if more than 75% of the field of view is covered)
                        mo.collided = 1;                                    % Mark that collision happened. 
                        mo.justCollided = 1;                                % Let the program outside know that it JUST happened
                    else
                        mo.justCollided = 0;
                    end
                    if(sumRet==0); mo.emptyFrames = mo.emptyFrames+1;       % Empty retina - maybe the object has left. Update the empty frames counter.
                    else; mo.emptyFrames = 0; end;                          % Not empty retina; reset the empty frames counter.
                else                                    % Testing stimuli, not training.
                    if(mo.d<0)                          % If reached full field...
                        mo.d = 0.00001;                 % freeze like that, in full field...
                        mo.dstep = 0;                   % and stop any further motion, not to have OFF signal from the retina
                    end
                end
                if(mo.shuffle); s = s(randperm(length(s(:)))); end;          % Randomize pixels if needs to be reshuffled                
                if(mo.noise>0); s = xor(s,floor(rand(size(s))+mo.noise)); end;
                % s = mo.retina; warning('Retina processing is turned off'); % Uncomment for visual testing reasons only, to see how the object behaves
            else
                s = zeros(dim);                                             % After collision - nothing to show
            end
    end
end
end


function h = myplot(varargin)
% h = myplot(Data)
% h = myplot(X, Y, Data) 
%
% Technical function.
% Draws a neat Heatmap, returns the handle for the surf object it created.

arguments = varargin(1:end);
switch length(arguments)
    case 1
        Data = arguments{1};   
    case 3
        x = arguments{1};
        y = arguments{2};
        Data = arguments{3};           
end        
h = surf((1:size(Data,2)+1)-.5, (1:size(Data,1)+1)-.5, zeros(size(Data,1)+1, size(Data,2)+1), Data);        

set(gca,'XLim',[0.5 size(Data,2)+0.5]);
set(gca,'YLim',[0.5 size(Data,1)+0.5]);
view(gca,2);
if(length(arguments)>=3)
    set(gca,'XTick',1:length(x),'XTickLabel',x);
    set(gca,'YTick',1:length(y),'YTickLabel',y);
end
colormap('hot');
useMap = get(gcf,'ColorMap');
colormap(gca,1-useMap);
shading(gca,'flat');

end


function M = remember(M,name,value)
% M = remember(M,name,value)
% Adds one more value to the list, and if it's the first row, also generates a title.

if(~isfield(M,'title')); M.title = {}; end;         % If no title field, create it
if(~isfield(M,'data'));  M.data = {}; end;         % If no data  field, create it
if(~isfield(M,'row'));   M.row = []; end;         % If no data  field, create it

if(~ismember(name,M.title))
    M.title{length(M.title)+1} = name;              % If it is a new title, add it
end
if(isempty(M.row))
    M.row = num2str(value);
else
    M.row = [M.row ',' num2str(value)];
end

% fprintf('%s \t %s\n',name,num2str(value));

end


function M = push(M)
% push(M)
% Saves this line, goes to the next line.

M.data{length(M.data)+1} = M.row;   % Add new row to the data
M.row = [];                         % Flush the row

end