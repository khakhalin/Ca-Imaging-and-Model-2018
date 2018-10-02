function varargout = model_stdp_multisens(stimType)
% model_stdp_multisens('full')
% [w,sel,th] = model_stdp_multisens_1()
%
% A simple model to play with plastic networks.
%
% Best way to generate lots of inputs, call it like that:
% for(i=1:50); model_stdp_multisens('looming'); end;

% Aug 08 2017: Toy model.
%   Transitioned to version _1 with inputs instead of forced random spiking.
% Aug 10 2017: Mostly functional
% Aug 14 2017: Blurry inputs, inactivation.
% Aug 15 2017: reporting, reshuffling, threshold adjustment stage
% Aug 17 2017: more varied collisions, more output measures
% Sep 28 2017: new outputs
% Jun 26 2018: Work on proper outputs started.
% Jul 03 2018: + looming and loomsh
% Sep 30 2018: General large-scale clean-up


%%% ------------------------    ------------ Constants ------------------------------------
% rng(3);                                       % Set random generator seed: keep commented for normal runs

saveResult = 1;                                 % If set to 1, runs a model and saves the results.
outFolder = 'C:\Users\Arseny\Documents\7_Ca imaging\Model outputs\';          % Where to save outputs

if(nargin<1); stimType = 'looming'; end;        % Stimulus type for training. Options: 
                                                % full      visual + mechanosensory hit-to-all at the end if it was a collision. Random mix of oblique transitions and collisions
                                                % vis       visual only, without a multisensory slam at the end
                                                % shuffle   like "vis", but scrambled. Note that sometimes hit-to-all at the end still happens because some vis stims have it
                                                % rand      fully random noise
                                                % looming   realistic looming stim (like vis, but only those that lead to a collision)
                                                % loomsh    exactly like "looming", but scrambled
                                                % crash     for testing: slightly randomized linear expanding (not realistic) central crashes
                                                % scramble  for testing: exactly like "crash", but scrambled
                                                % trans     just a translation. Currently only half-implemented: not randomized
if(nargin<1); showFigures = 1; else; showFigures = 0; end;  % Don't generate figures if called externally (saves time)

brainDim = 9;                                   % Size of the square matrix of cells (and retina). 9 or higher is recommended
nCells = brainDim^2;                            % Number of neurons (any number from 5 about 50. Above that it needs more noise, or more time to converge)
nStim = 500;                                    % Approximate number of stimuli (say 2000)
nTick = 25;                                     % Number of time ticks in each stimulus (100 for 1s, but 25 is good). Should be > than collision time
tmax = nStim*nTick;                             % max time the model will actually run (should be about 1k-5k); runs with 10 ms steps
nTickRecalibrate = 2000;                        % Number of ticks (not stim) for network overall excitability (th) recalibration before testing. 1000?

% --------- sTarget : target avereage spiking rates for each cell. Several options here:
% sTarget = ones(nCells,1);                     % Identical. 
sTarget = randn(nCells,1)/5+1;                  % Noisy
% sTarget = 1./(1+rand(nCells,1)*5);            % Realistic: Most cells are non-spiky. Param of 5 seems good.
sTarget = sTarget/mean(sTarget)*5/nCells;       % Normalization. If sTarget==k/nCells, on average k neurons will be active at any moment (here k=5)
% figure; histogram(sTarget,20); return;         % Simple test to see what we got

sAvDecay = 0.95;                                % Decay coeff for running averaging of spiking = 1-tstep/tau; for tau=200ms k=0.95 (seems to be good)
thSens = 0.1;                                   % Threshold change rate in the beginning. 0.1 is good
plasticityType = 'STDP';                        % Two options: STDP and Hebb
kHebb = 0.25;                                   % STDP change rate. 1 is too fast (for nstim=1000), 0.1 takes too much time. Optimal: 0.25
blurDistance = 0.0;                             % Blurring distance for RGC inputs. Set to 0 for no blur, 0.5 for weak, 1 for pretty strong
incompleteProjections = 0.0;                    % Amount of noise to add to retinal inputs (0 for none, 1 for full rand)
inactivationMode = 'none';                      % Options: 'none' (recommended), 'fast' (10 ms period), and 'slow' (~200 ms, based on sAvDecay)
competitionMode = 'slide';                      % How synaptic competition os realized. Options: in, out, global, decay, none, softin, softout, satin, satout, slide (best)
totalWeight = 1.5;                              % The forced sum of all outputs (aka How many saturated outputs could a cell have). Default: 1.5
synapticDecay = 1-0.001;                        % Synaptic decay constant; is only relevant if decay is used for synaptic competition. Something like 1-0.001 may work

% Set that kinda works (2017): brainDim  9, nStim 1000, nTick 25, decay 0.98, thSens 0.1, kHebb 0.2, no blur (0), no inactivation
% Another set (2018): brainDim 9, nStim 500, nTick 25, decay 0.95, thSens 0.1, kHebb, 0.25, noblur, no inactivation
% Random stimulation ('rand') produces slowest convergence (stabilization) of network properties, while looming ('crash') are the fastest.
% So when optimzing a set of parameters (especially nStim and kHebb) make sure it kind of works for both

fprintf('Starting a model:%d cells, %d stimuli of %s type, %s rule, %s competition, %s inactivation\n',...
    nCells, nStim, stimType, plasticityType, competitionMode, inactivationMode);

%%% ------------------------------------ Starting parameters ------------------------------------
w = rand(nCells);                               % Random connections
% w = rand(n).*floor(rand(n)+0.95);             % Shotgun random connections if necessary
w = w-diag(diag(w));                            % remove self-connections
%w = bsxfun(@times,w,1./sum(w,2));               % Scale the sum of outputs for every neuron to one
% w = zeros(n); for(q=2:n); w(q-1,q) = 1; end; w(n,1) = 1; % This row creates one giant cycle, for testing purposes
th = 1./sTarget/nCells + 0.1*rand(nCells,1);    % A noisified guess about good original thresholds
[meshx,meshy] = meshgrid(1:brainDim,1:brainDim);          % Cell positions, for figures and calculations
hiddenWeight = 0.5*ones(nCells,1);              % A vector of hidden links to (or from) cells; used in soft types of normalization

%%% ------------------------------------ Technical variables ------------------------------------
wstart = w;                                     % keep starting weights as a memento, to plot at the end
s = floor(rand(nCells,1)+sTarget);              % spikes vector (start with something)
s1 = zeros(nCells,1);                           % previous spikes (none)
sAv = sTarget;                                  % recent average spiking. Assume all spiked ideally at first
sHist = zeros(nCells,tmax);                     % history of spikes
sAvHist = zeros(nCells,tmax);                   % history of averaged spikes
thHist = zeros(nCells,tmax);                    % history of thresholds
wHist = zeros(nCells,tmax);                     % history of synaptic inputs to one neuron
slider = zeros(tmax,1);                         % Easing coefficient
U = [];                                         % Future output structure, to be saved

blur = eye(nCells);                             % Blurring retinal inputs
shotgun = ones(nCells)*(1-incompleteProjections) + rand(nCells)*incompleteProjections;
if(blurDistance>0)
    for(i=1:nCells); for(j=1:nCells); blur(i,j) = exp(-((mod(i-1,brainDim) - mod(j-1,brainDim))^2 + (floor((i-1)/brainDim)-floor((j-1)/brainDim))^2)/(2*blurDistance^2)); end; end;
    blur = blur.*shotgun;
    blur = bsxfun(@times,blur,1./sum(blur,2));  % Practical normalization - a bit weird at the edges
    blur = bsxfun(@times,blur,1./sum(blur,1));
else
    blur = blur.*shotgun;
end

%%% ------------------------------------ Technical testing area ------------------------------------
% Uncomment tests here for bypassing the main loop and running these tests instead.

% figure; subplot(1,2,1); myplot(reshape(blur(floor(nCells/3),:),dim,dim)); subplot(1,2,2); plot(blur(floor(nCells/3),:),'.-'); return;
% test_inputs(brainDim,stimType,blur); return;         % Draws a few stimuli, as they would be presented to the network.                                                 

%%% ------------------------------------ Main loop ------------------------------------
t = 0;
si = 0;
nCollisions = 0;
iStage = 0;                                     % Current developmental stage
counter = 0;                                    % Counter for periodical network assessment (at end of this loop)
dataCollected = [];                             % Here data will be collected
timerOriginal = tic;
timer = tic;
while(t<=tmax)                                  % We control for the total number of ticks, and not for the number of stimuli
    if(toc(timer)>10)                           % Every 10 seconds, update the console (in case the simulation gets real slow)
        timer = tic;
        if(showFigures); fprintf('...now at tick=%d, stimulus=%d (about %d seconds to go)\n',t,si,floor(toc(timerOriginal)/t*(tmax-t))); end
    end
    si = si+1;                                  % New stimulus
    [mo,~] = generateInput(brainDim,stimType);  % Create motion object for this stimulus
    for(ti=1:nTick)                             % Run this stimulus
        t = t+1;                                % Time step
        slider(t) = exp(-7*t/tmax);             % Gradually decreasing value, for easing (trasitioning) functions. 5 is fine
        s1 = s;                                 % History of spiking (spikes from the previous run)
        switch inactivationMode
            case 'none'; s = w*s;                                	% Excitatory transmission, no refractory period
            case 'fast'; s = w*s./(1+s);                            % Excitatory transmission + short inactivation period (10 ms)
            case 'slow'; s = w*s./(1+sAv);                          % Excitatory transmission + long inactivation period (~200 ms)
        end
        
        [mo,input] = generateInput(mo);                             % Sensory input
        s = s + blur*input(:);
        if(mo.justCollided)
            if(strcmp(stimType,'full')); input = ones(size(input)); end;    % Slam with multisensory input
            nCollisions = nCollisions+1;                            % Count collisions. ('Stimulus %d collided at tick %d\n',si,ti); 
        end        

        % --- Facilitate transition and avoid bottlenecks:
        % s = s*(1*(1-slider) + 3/sum(s)*slider(t));    % Forceful inhibition. Doesn't work that well.
        % s = s/(1+max(0,sum(s)-dim));                  % A better attempt at forced inhibition. Kicks in when more than dim cells are active
        s = s/(1+max(0,sum(s)-brainDim)*slider(t));     % Same as above, but does not stay forever; instead it gradually disappears (because of "slider")
        %s = s+randn(ncells,1)/3*slider(t);             % Gaussian noise to shake the system early on

        s = logisticf(s,th);                    % Spiking (soft non-linear logistic function, defined below)    
        sAv = sAv*sAvDecay + s*(1-sAvDecay);    % Calculate running average spiking
        th = th - (sTarget-sAv)*thSens;         % Homeostatic plasticity that stays the same

        %%% Synaptic plasticity
        switch plasticityType
            case 'STDP'
                %w = w+diag(s-s1)*w*diag(s1)*kHebb;         % STDP: reward if s1 interacted (w) with s, but a punishment if with s1. This formula is very slow though        
                w = w+kHebb*bsxfun(@times,bsxfun(@times,s(:)-s1(:),w),s1(:)'); % STDP; Exactly the same formula as above, but ~100 times faster (literally)
            case 'Hebb'
                w = w+kHebb*bsxfun(@times,bsxfun(@times,s(:),w),s1(:)'); % Simple Hebb (as an alternative)
            otherwise
                error('Wrong plasticity type');
        end
                
        %%% Synaptic competition (several alternatives; pick any one)
        switch(competitionMode)
            case 'out';     w = bsxfun(@times,w,totalWeight./sum(w,1));    % Aggressive synaptic scaling of all OUTputs of every neuron (hard renormalization at every step)
            case 'in';      w = bsxfun(@times,w,totalWeight./sum(w,2));    % Aggressive synaptic scaling of all INputs of every neuron (hard renormalization at every step)
            case 'softout'; temp = (sum(w,1)+hiddenWeight(:)'); w = bsxfun(@times,w,totalWeight./temp); hiddenWeight = totalWeight*hiddenWeight./temp'; % Careful with dimensions
            case 'softin';  temp = (sum(w,2)+hiddenWeight(:));  w = bsxfun(@times,w,totalWeight./temp); hiddenWeight = totalWeight*hiddenWeight./temp;  % In this mode cells have "hidden synapses" that compete with visible ones
            case 'decay';   w = w*synapticDecay;                        % Simple synaptic decay: very sensitive to the decay value; not pleasant
            case 'global';  w = w/sum(w(:))*nCells*totalWeight;         % Limiting total number of edges to the number of cells (because every w_ij is <=1), adjusted by totalWeight
            case 'satin'                
                w = bsxfun(@times,w,1./max(1,max(w,[],2)));             % Normalizing only if weights are creeping above 1
                %w = bsxfun(@times,w,sum(w>1,2)
            case 'satout'
                w = bsxfun(@times,w,1./max(1,max(w,[],1)));             % Normalizing only if weights are creeping above 1
                %w = bsxfun(@times,w,sum(w>1,2)
            case 'slide';                                               % Gradual synaptic scaling. Danger of oscillations, needs sliding parameters
                temp1 = bsxfun(@times,w,totalWeight./sum(w,1));         % trimmed outs
                temp2 = bsxfun(@times,w,totalWeight./sum(w,2));         % trimmed ins
                w = 0.4*w + 0.3*temp1 + 0.3*temp2;                    % The sliding equation - if too slow, hebb may get out of control still.
            case 'none';    % nothing
            otherwise; error('Illegal synaptic competition mode.');
        end
        
        w = max(0,min(1,w));                        % Clip negative weights, and those that are above 1
        w = w-diag(diag(w));                        % Keep self-connections at 0
        
        thHist(:,t) = th;                       % Store threshold history for plotting
        wHist(:,t) = w(:,1);                    % Weight history for one neuron (for troubleshooting)
        sHist(:,t) = s;                         % History of spiking
        sAvHist(:,t) = sAv;                     % History of "recent average spiking"
        
        if(t==1 || mod(t,1000)==0)              % Every 1000th time step - calculate some stats and save them
            counter = counter+1;
            newDataLine = [];            
            temp = efficiency_wei(w');                      newDataLine = [newDataLine temp];
            temp = assortativity_wei(w',4);                 newDataLine = [newDataLine temp]; % 1 out-in, 2 in-out, 3 out-out, 4 in-in
            [~,temp] = modularity_dir(w');                  newDataLine = [newDataLine temp];
            temp = mean(clustering_coef_wd(w'));            newDataLine = [newDataLine temp];
            [flow, revFlow] = network_flow(w');             newDataLine = [newDataLine flow revFlow];  % Global reach but via Katz centrality
            %cycl = myCyclicity(w');                         newDataLine = [newDataLine cycl]; % It drops so fast that there's no point to build it like that
            dataCollected = [dataCollected; t/1000 newDataLine];
            %dispf(newDataLine,'%10.5f');
        end
        if(saveResult && (t==10 || t==floor(tmax/4) || t==floor(tmax/2) || t==floor(tmax*3/4) || t==tmax))  % Save current developmental stage
            iStage = iStage+1;
            fprintf('Got to stage %d (%d%%)\n',iStage,round(t/tmax*100));
            U.stage(iStage).w = w;        % Matrix of weights
        end
        
        if(mo.justCollided); break; end;        % After a collision the rest of the trace will be empty for sure, so no need to wait. Start a new one.
        if(mo.emptyFrames>3); break; end;       % The object has left the screen, no need to wait. Start a new one.
    end
end
dataLabels = {'Efficiency','InInAssort','Modular','Cluster','Flow','revFlow'};       % The order in which dataCollected stores measurements
nStim = si;                                     % Actual number of stimuli that happened (in case si is overwritten later)

if(saveResult)                                      % If it's time to save the results
    U.measurements = dataCollected;                 % Just save all measurements bulk
    U.type = stimType;                              % Now save all key parameters that the analysis tool needs to know
    U.sTarget = sTarget;                            % This one we'll need for testing
    U.blur = blur;                                  % RGC->OT connections matrix; is also needed for testing
    U.totalOut = totalWeight;                       % This one and below we'll need for accounting purposes (to know what is what)
    U.sAvDecay = sAvDecay;
    U.thSens = thSens;
    U.plasticityType = plasticityType;              % STDP or simple Hebb
    U.kHebb = kHebb;
    U.blurDistance = blurDistance;
    U.incompleteProjections = incompleteProjections;
    U.inactivationMode = inactivationMode;          % Remember options
    U.competitionMode = competitionMode;
    
    fileName = [datestr(now,'yymmddHHMMss') '-' competitionMode '-' stimType];
    save([outFolder fileName '.mat'],'U');
end

if(showFigures) % --- How everything changes (or not) with time (which is stored in dataCollected(:,1))
    figure('Color','white'); 
    for(q=1:length(dataLabels))
        subplot(2,4,q); plot(dataCollected(:,1),dataCollected(:,1+q),'.-'); title(dataLabels{q});
    end
end


%%% ------------------------------------ Figures ------------------------------------
if(showFigures) % --- Main figure (mugshot)
    figure('Color','white');
    subplot(2,3,1); set(gca,'FontSize',8); myplot(wstart(end:-1:1,:)); title('Starting weights');
    subplot(2,3,2); set(gca,'FontSize',8); myplot(     w(end:-1:1,:)); title('Final weights');
    %subplot(2,3,3); gplot(round(w),[cos((1:ncells)'/ncells*2*pi) sin((1:ncells)'/ncells*2*pi)],'.-'); axis('off'); title('Graph');    
    subplot(2,3,3); gplot(round(w),[meshx(:) meshy(:)],'.-'); axis('off'); title('Graph');
    subplot(2,3,4); plot(1:size(wHist,2),wHist'); set(gca,'FontSize',8,'LineWidth', 1.2); title('w_{i1} history');
    subplot(2,3,5); plot(1:size(thHist,2),thHist'); set(gca,'FontSize',8,'LineWidth', 1.2); title('Threshold history'); hold on; plot(1:length(slider),slider,'r.-'); hold off;
    subplot(2,3,6); myplot(resample(sHist',1,ceil(size(wHist,2)/50))'); hold on; plot(resample(sum(sHist),1,ceil(size(wHist,2)/50)),'k-','LineWidth',2); hold off; title('Spike history');
    %myplot(sAvHist); title('Av Spike history'); % plot(sAvHist'); 
    drawnow();
end

return; % --- Uncomment if no network testing is needed





%%% --------- Recalibration of thresholds before actual testing ---------
% Helps to get rid of short-term effects of most recent learning. Like sleeping overnight. 
% Makes it fair, levels the field for every model, especially after training on different stimuli.
% (Sometimes the most recent stimulus was intense and all neurons are a bit instinsically suppressed, which would bias the test).
% But no hebbian plasticity here; only the thresholds, recalibrated on noise.
for(q=1:nTickRecalibrate)                                       % Recalibrate thresholds only, leaving w intact
    switch inactivationMode
        case 'none'; s = w*s;                                	% Excitatory transmission, no refractory period
        case 'fast'; s = w*s./(1+s);                            % Excitatory transmission + short inactivation period (10 ms)
        case 'slow'; s = w*s./(1+sAv);                          % Excitatory transmission + long inactivation period (~200 ms)
    end
    input = floor(rand(size(s))+1/brainDim);                         % Keep the noise here rather weak (dim/ncells) to enable chain-ringing
    s = s + blur*input(:);
    s = logisticf(s,th);                    % Spiking (soft non-linear logistic function, defined below)    
    sAv = sAv*sAvDecay + s*(1-sAvDecay);    % Calculate running average spiking
    th = th - (sTarget-sAv)*thSens;         % Homeostatic plasticity that stays the same
end


%%% ------------------------------------ Testing selectivity ------------------------------------
typeCycle = {'flash','scramble','crash'};
s = sTarget;                                    % spikes vector set to ideal (needs to be reset after training)
sAv = sTarget;                                  % recent average spiking. Assume all spiked ideally at first
nTestReps = 10;
nStimTypes = 3;
for(si=1:nStimTypes)                                     % For each testing stimulus. 1=F, 2=C, 3=transition
    spikeHistory{si} = zeros(nCells,nTick);              % To collect spiking
    inputTrace{si} = zeros(1,nTick);
    spikeOutput{si} = [];
    for(q=1:nTestReps)                          % Repeat several times
        [mo,~] = generateInput(brainDim,typeCycle{si});        
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
            spikeHistory{si}(:,ti) = spikeHistory{si}(:,ti)+s;          % Remember spiking by each neuron at every moment of time
        end        
        spikeOutput{si} = [spikeOutput{si} sum(spikeHistory{si},2)];    % Total spiking for this neuron
    end
    spikeHistory{si} = spikeHistory{si}/nTestReps;    
    inputTrace{si}   = inputTrace{si}/nTestReps;
end

resF = mean(spikeOutput{1},2);                          % Average total output for flashes. A column of numbers
resS = mean(spikeOutput{2},2);                          % Scramble
resC = mean(spikeOutput{3},2);                          % Crash
% sel = (resC-resF)./resF;                              % A simple selectivity measure (0 for no selectivity)
selFC = (resC-resF)./std([bsxfun(@plus,spikeOutput{1},-resF) bsxfun(@plus,spikeOutput{2},-resC)],[],2); % Cohen's d

if(showFigures) % --- How selectivity interacts with cell place within the network
    figure('Color','white');
    lineStyle = {'k-','r-','b-'};
    %h6 = subplot(3,3,6); hold on;
    for(si=1:length(spikeHistory))
        subplot(3,3,si);   myplot(reshape(mean(spikeOutput{si},2),brainDim,brainDim)); title(typeCycle{si}); colorbar(); %caxis([0 2]);
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
    distanceToCenter = sqrt((meshx(:)-brainDim/2).^2 + (meshy(:)-brainDim/2).^2);
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



function test_inputs(dim,type,blur)
nrows = 18;
ncols = 25;
figure;
for(si=1:nrows)
    [mo,~] = generateInput(dim,type);        
    for(ti=1:ncols)
        [mo,s] = generateInput(mo);
        input = blur*s(:);
        %if(mo.justCollided); input = ones(size(s)); end;
        axes('units','norm','pos',[(ti-1)/ncols (si-1)/nrows 0.9/ncols 0.9/nrows]); % Subplot is too slow here
        set(gca,'XtickLabel','','YTickLabel','','nextplot','add');
        myplot([],[],reshape(input,dim,dim));        
        caxis([0 1]);
        if(mo.justCollided); break; end;    % Done, via collision        
        if(mo.emptyFrames>3); break; end;   % Done, left the field
    end
end
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
            mo.dfinal = -3;                                     % Final distance (doesn't make much sense in this case). Adjusted to make the full transition exactly 2 frames
            mo.startAngle = 0;                                  % Starting offset angle, from the center of view
            mo.moveAngle = 0;                                   % Moving direction (opposite to the offset angle + some noise)
            mo.start = [dim/2, dim/2]+0.5;                      % Starting position (doesn't matter in this case?)
            mo.v = 0;                                           % Lateral speed
            mo.shuffle = 1;                                     % This one needs to be shuffled, to have some minimal variability in responses
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
            mo.moveAngle = mo.startAngle+pi+(rand(1)-0.5)*pi/3; % Moving direction (opposite to the starting angle + noise)
            offset = rand(1)*(dim/2+mo.r-1);                    % Starting offset from the center (from 0 to not visible)
            mo.start = [dim/2 + cos(mo.startAngle)*offset, dim/2 + sin(mo.startAngle)*offset]+0.5;  % Starting position
            mo.v = dim/mo.period;                               % Speed. It will cross the field in this time
            if(strcmp(mo.type,'shuffle'))
                mo.shuffle = 1;                                 % Set Shuffle to 1 if necessary
            end
        case {'looming','loomsh'} % Training transitions that always end with a crash            
            mo.dstart = 1;                                      % Starting distance to the object (the same as all)
            mo.dfinal = 0.01;                                   % Final distance to the object (almost 0, but not quite)
            mo.startAngle = rand(1)*2*pi;                       % Starting offset directoin relative to the center of the screen
            offset = rand(1)*(dim/2+mo.r-1);                    % Starting offset from the center (from 0 to not visible)
            mo.moveAngle = mo.startAngle+pi+(rand(1)*2-1)*atan(offset/mo.r); % Moving direction (opposite to the offset direction + some leeway)
            mo.start = [dim/2 + cos(mo.startAngle)*offset, dim/2 + sin(mo.startAngle)*offset];  % Starting position
            mo.v = (offset+(rand(1)*2-1)*mo.r)/mo.period;       % Lateral speed; made lower than for "vis" group, to keep collisions more central
            if(strcmp(mo.type,'loomsh'))
                mo.shuffle = 1;                                 % Set Shuffle to 1 if needed
            end
            mo.period = mo.period + floor(rand(1)*5);           % Make some stimuli a bit slower (to have them a bit more different)
        case 'rand' % Just fill with random noise
            mo.dstart = 1; mo.dfinal = 1; mo.startAngle = 0; 
            mo.moveAngle = 0; mo.start = [0 0]; mo.v = 0; 
        otherwise
            error(sprintf('Unrecognized stimulus type: %s',mo.type));
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
    s = zeros(dim);                                     % Even though it's not used at init, populate with something to prevent from crashing
    
else  % --- Generate new frame using motion object "mo"
    mo = varargin{1};
    ncells = numel(mo.retina);
    dim = sqrt(ncells);    
    
    switch mo.type
        case 'rand'
            noiseLevel = 1/dim;      % Noise level.
            s = floor(rand(ncells,1)+noiseLevel);
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
                mo.rgc = xor(mo.retina,oldretina);                          % Only report change (on/off cells) - proper mode, as it should be
                % mo.rgc = mo.retina;                                       % Image itself: is only good for stimulus debugging; wouldn't work in a real model!
                % mo.rgc = mo.rgc*rgcDecay + xor(mo.retina,oldretina);      % Older version with decay; doesn't work (see comments)                
                s = mo.rgc;
                
                if(strcmp(mo.type,'full') || strcmp(mo.type,'vis') || strcmp(mo.type,'shuffle') ...
                        || strcmp(mo.type,'looming') || strcmp(mo.type,'loomish')) % One of the training stimuli: detect if collision had happened
                    sumRet = sum(mo.retina(:));
                    if(sumRet>0.75*ncells & mo.collided==0)                 % Collision happened (if more than 75% of the field of view is covered)
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
                if(mo.shuffle); s = s(randperm(length(s(:)))); end;          % Randomize pixels if needs to be reshuffled.              
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
