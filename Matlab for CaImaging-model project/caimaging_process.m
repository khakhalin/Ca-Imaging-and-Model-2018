function res = caimaging_process(S)
% CAIMAGING_PROCESS(S)

% Dec 17 2013: Created.
% Dec 18 2013: Improved.
% Dec 30 2013: + more visualizations.
% Jan 31 2014: + filtered averaged output. + More careful approach to the time scale.
% Mar 04 2014: + comparison of cells across stimuli.
% Mar 05 2014: + Comparison of latencies (weighted over time) and selectivities.
% Mar 06 2014: Tiding up.
% May 01 2014: + History of selected cells.

%%% --------------------------------------------------------------- HEADER ---------------------------------

doAverage = 1;          % Average response, over time
doAverageMap = 0;       % Spatial map of response strength
doLatencyMap = 0;       % Latency map
doHeatMaps = 0;         % History of total responses by cell
doHistoryOf = [];       % If not empty, for these cells some more detailed history will be given.
doCompareCells = 0;     % Compare cells in their responses to different stimuli
doPca = 0   ;           % PCA

        % Below: Latency adjustments, in ms. These supposedly should come from the videos, but this approach 
        % totally don't work, as the eye is (apparently) more sensitive than the video-based curve.
        % Also it is not necessarily helpful when "Crash" response starts so late compared to "Flash" 
        % (unless you are trying to assess the actual delay in the retina).

%%% % key = 'fcgb'; latAdj = [307 490 490 775];     
% key = 'ofsc'; latAdj = [490 307 490 490];           % Visual inputs CALCULATED. latAdj goes in ms.
% key = 'ofsc'; latAdj = [800 360 800 680];           % Visual inputs EMPIRICAL. latAdj goes in ms.
key = 'cfs'; latAdj = [500 300 700];           % EMPIRICAL
% latAdj = [1 1 1 1]*0;             % For debugging
% key = 't'; latAdj = [290];        % Tentacle. I don't really know the delay yet, so it's guessing.

% myColorMap = [linspace(245,173,50) , linspace(173,49  ,50) ; ...
%               linspace(252,221,50) , linspace(221,163  ,50) ; ...
%               linspace(185,142,50) , linspace(142,84,50)]'/256;

myColorMap = [linspace(247,120,50) , linspace(120,0  ,50) ; ...
              linspace(252,198,50) , linspace(198,104  ,50) ; ...
              linspace(185,121,50) , linspace(121,55,50)]'/256;
                                    
nTrials = length(S);
[n,nCells] = size(S(1).dataS);                      % n = ntime
for(iTrial=1:nTrials)
    n = max(n,size(S(iTrial).dataS,1));             % In case some traces are 1-2 points longer than the others due to sketchy smoothing of gaps (by caimaging_basic.m)
end
nType = length(key);
dt = (S(1).timeS(end)-S(1).timeS(1))/n;             % Time step in seconds
time = (1:n)*dt*1000;                               % in ms
liveZone = 2;                                       % 2 seconds
ntimeticks = round(liveZone/dt);                    % N points that make the liveZone. (most probably mean(diff(time))==time(2)-time(1), but who knows...
fType = find(key=='f');                             % To be able to compare 'flash' to 'crash' in different protocols.
cType = find(key=='c');                             %   BEWARE, there's no workaround for non-fc-containing data (yet).
if(isempty(fType*cType))
    fprintf('Either F or C stimuli are absent form the key. Prepare for errors below!\n');
end

x = S(1).xy(:,1);
y = S(2).xy(:,2);  
y = max(y)-y;           % Inversed, as Y is measured right-down by NES instruments, while I like it right-up

res = [];               % In case any results are brought to the console



%%% --------------------------------------------------------------- PROCESSING ---------------------------------

%%% ---------- History of individual cells
relativeThreshold = 0.99;   % Threshold for spiking
if(~isempty(doHistoryOf))
    for(iCell=doHistoryOf)
        figure;
        subplot(1,2,1);
        bag = [];        
        for(iTrial=1:nTrials)
            bag = concatenan(bag,S(iTrial).dataF(:,iCell));
        end
        bag = bsxfun(@plus,bag,-bag(1,:));
        myplot(bag);  title(sprintf('Cell %d: F-data',iCell));
        subplot(1,2,2);    
        if(1)
            [b,a] = butter(3,1/5);
            plot(filter(b,a,bag(:,[1 4 7 37]))); title('F-data, as plots');
        else
            bag = [];
            for(iTrial=1:nTrials)
                bag = concatenan(bag,S(iTrial).dataS(:,iCell));
            end
            spikeThreshold = quantile(bag(:),relativeThreshold);
            myplot(bag.*(bag>spikeThreshold));  title('S-data, relative threshold');                   
        end
    end
end


%%% ---------- Simple unthresholded averages, over time (AV), and total (AMPSRAW)
for(iType = 1:nType)
    avRaw{iType} = zeros(n,nCells);
    typeSweepsCount(iType) = 0;
end
ampsRaw = nan(nTrials,nCells);
for(iTrial = 1:nTrials)    
    iType = mod(iTrial-1,nType)+1;
    temp = S(iTrial).dataS;
    if(size(temp,1)<n)      % Unify length
        temp = [temp; zeros(n-size(temp,1),nCells)];
    end
    ampsRaw(iTrial,:) = sum(temp,1);                            % Raw amplitudes (noisy, bleaching, terrible)
    avRaw{iType} = avRaw{iType} + temp;
    typeSweepsCount(iType) = typeSweepsCount(iType)+1;
end
%%% Pre-processing: cut the edges and synchronize time-scales
for(iType=1:nType)
    ind = find(time>=latAdj(iType),1,'first');
    avRaw{iType} = avRaw{iType}(ind-1 + (1:ntimeticks),:);            % A weird way to cut, used because time scales are not quite uniform    
    avRaw{iType} = avRaw{iType}/typeSweepsCount(iType);    
end
time = time(1:size(avRaw{1},1))';


%%% --------------- Thresholded averages and total amplitudes of response (AMPS and its friends)
%%% ---- Approach #1 : per-trial thresholding (ampsPerTrial)
relativeThreshold = 0.99;                                           % Quantile used to detect spikes (0.98 works fine; 0.99 may be a bit darkish)
for(iTrial = 1:nTrials)
    iType = mod(iTrial-1,nType)+1;                                  % Stimulus type for this trial, for local use
    type(iTrial) = iType;                                           % Stimulus type for this trial, for later use
    indFirst = find(time>=latAdj(iType),1,'first');                 % Adjust the latency    
    tempTraces = S(iTrial).dataS(indFirst-1 + (1:ntimeticks),:);    % Temporary variable with all traces for this trial
    spikeThreshold = quantile(tempTraces(:),relativeThreshold);     % Global spike threshold (for all cells simultaneously)
    tempThresholded = tempTraces;                                   % Create a copy.
    tempThresholded(tempThresholded>=spikeThreshold) = 1;           % This is questionable actually. Can we have more than one spike within same 10 ms bin?
    tempThresholded(tempThresholded< spikeThreshold) = 0;    
    ampsPerTrial(iTrial,:) = sum(tempThresholded,1);                % Globally thresholded amplitudes.    
end
clear tempTraces tempThresholded spikeThreshold;                    % Clear temp values to avoid confusion


%%% ---- Approach #2 : per-cell thresholding (ampsPerCell)
relativeThreshold = 0.99;                                           % Now prepare per-cell thresholded data (ampsPerCell). Lower threshold (otherwise F isn't visible)
for(iCell=1:nCells)                                                 
    for(iTrial = 1:nTrials)
        tempTraces(:,iTrial) = S(iTrial).dataS(indFirst-1 + (1:ntimeticks),iCell);  % Collect all data for this cell across all trials        
    end
    spikeThreshold = quantile(tempTraces(:),relativeThreshold);
    for(iTrial=1:nTrials)
        ampsPerCell(iTrial,iCell) = sum(tempTraces(:,iTrial)>=spikeThreshold);
    end
end
clear tempTraces;


%%% ---- Approach #3 : per-group of trials thresholding (amps)
for(iType=1:nType)
    av{iType} = zeros(ntimeticks,nCells);                               % Placeholder for average traces
end
relativeThreshold = 0.99;                                               % Quantile used to detect spikes (0.98 works fine; 0.99 may be a bit darkish)
amps = nan(nTrials,nCells);                                             % Here per-trial thresholded amplitudes will be stored
for(iTrial = 1:nTrials)
    iType = type(iTrial);
    indFirst = find(time>=latAdj(iType),1,'first');                     % Adjust the latency    
    tempTraces{iTrial} = S(iTrial).dataS(indFirst-1 + (1:ntimeticks),:);% Temporary variable with all traces for this trial
    if(mod(iTrial,nType)==0)                                            % Full cycle of stimuli is over    
        bag = [];
        for(iAdd=1:nType)
            bag = [bag; tempTraces{iTrial-(iAdd-1)}(:)];                % Collect all data here
        end
        spikeThreshold = quantile(bag,relativeThreshold);               % Global-ish spike threshold (for all cells simultaneously, within these several trials)
        for(iAdd=1:nType)
            tempThresholded = tempTraces{iTrial-(iAdd-1)};              % Create a copy of traces.
            tempThresholded(tempThresholded>=spikeThreshold) = 1;       % Limit spikes
            tempThresholded(tempThresholded< spikeThreshold) = 0;
            tempTraces{iTrial-(iAdd-1)} = tempThresholded;              % Store for the future
            av{type(iTrial-(iAdd-1))} = av{type(iTrial-(iAdd-1))}+tempThresholded;            
            amps(iTrial-(iAdd-1),:) = sum(tempThresholded,1);
        end        
    end
end
clear tempTraces tempThresholded spikeThreshold;                    % Clear temp values to avoid confusion
for(iType=1:nType)
    ampsByType(iType,:) = mean(amps(type==iType,:),1);
    av{iType} = av{iType}/sum(type==iType);
end

[~,iSortType] = sort(type);                                 % To sort responses by stimulus
% iSortType=1:length(iSortType);
[~,iSortAmps] = sort(sum(amps,1));                          % To sort cells by response amplitude
if(doHeatMaps)
    figure;
    myplot(ampsRaw(iSortType,iSortAmps));   title('Raw amplitudes');
    figure;                                     % Heatplot of cell responses in different trials. Global thresholding.
    myplot(amps(iSortType,iSortAmps)); title('Thresholding per group of trials');    
    figure;
    myplot(ampsPerCell(iSortType,iSortAmps)); title('Thresholding per cell');
end



centerCell = 53;    % Cell that lies about in the center of where C emerges and B ends. 33 for dataset from 131213
circR = 30;         % Radius of the circle to define "central cells"
xc = S(1).xy(centerCell,1);
yc = S(1).xy(centerCell,2);
indIn  = sqrt((x-xc).^2+(y-yc).^2) <= circR;



%%% ------------ Average response as f(t)
if(doAverage)
    timeSlice = 1:find(time>1000,1,'first');    % First second.  No need to use latAdj, as av{} are trimmed already.
    figure;
    for(iType = 1:nType)    
        subplot(nType,1,iType);
        hold on;
        plot(time(timeSlice),av{iType}(timeSlice,:),'-','Color',[7 7 8]/8);
        plot(time(timeSlice),mean(av{iType}(timeSlice,:),2),'b-');
        hold off;
        xlim([min(time(timeSlice)) max(time(timeSlice))]);
        ylim([0 0.6]);
        title(key(iType));
        drawnow;
    end
    
    figure;
    for(iType = 1:nType)    
        subplot(nType,1,iType);
        hold on;
        plot(time(timeSlice),avRaw{iType}(timeSlice,:),'-','Color',[7 7 8]/8);
        plot(time(timeSlice),mean(avRaw{iType}(timeSlice,:),2),'b-');
        hold off;
        xlim([min(time(timeSlice)) max(time(timeSlice))]);
        ylim([0 0.6]);
        title(key(iType));
        drawnow;
    end
end


%%% ---------- Compare cells selectivity
for(iCell=1:nCells)
   pSelective(iCell) = anova1(amps(:,iCell),type,'off');                                    % Overall selectivity across all stimuli
   [~,pSelectiveFC(iCell)] = ttest2(amps(type==fType,iCell),amps(type==cType,iCell));       % T-test between responses to C and F. 
                                                            % IT NOW USES AMPS = NORMALIZED PER TRIAL = MOST SELECTIVITY IS GONE!!!
   kSelectiveFC(iCell) = -log(pSelectiveFC(iCell))*sign(ampsByType(cType,iCell)-ampsByType(fType,iCell));  % Signed selectivity index. Positive if C>F, negative if not, ~zero is non-significant.
end
reallySelective = 10;                                       % Just a random "big number" to put instead of infinity for pval==1
kSelectiveFC(isinf(kSelectiveFC)) = reallySelective*sign(kSelectiveFC(isinf(kSelectiveFC)));
kSelectiveFC(isnan(kSelectiveFC)) = 0;                      % For some reason ttest2 returns NaN for identical inputs, intead of a 1 as I'd expected. and -log(1) = 0.
selectiveCells = pSelective < 0.05;                         % A list of selective cells.

if(doCompareCells)
    figure;                                     % Selectivity vs amplitude ------- CHANGE THIS TO RAW AMPLITUDE!!!!
    hold on; plot(sum(amps),kSelectiveFC,'bo'); 
    plot(get(gca,'XLim'),[1 1]*-log(0.05),'r--'); 
    [rho,pval] = corr(sum(amps)',kSelectiveFC');
    p = polyfit(sum(amps)',kSelectiveFC',1);        
    plot(get(gca,'XLim'),polyval(p,get(gca,'XLim')),'g-');
    hold off; 
    xlabel('Amplitude'); ylabel('Selectivity for F-C pair')   
    title(sprintf('Selectivity vs. amplitude. Correlation: r = %s, p = %s',myst(rho),myst(pval)));
    
    spearman_plot(sum(amps),kSelectiveFC);    
end


%%% ------------ Map of Average reponse amplitudes
if(doAverageMap)
    loC = min(mean(av{1},1));   % For color scale
    hiC = max(mean(av{1},1));
    figure;
    for(iType = 1:nType)
        if(nType>1)
            subplot(2,2,iType);
            title(key(iType));
            hold on;
        end        
        % scatter(x,y,15,colorscale(mean(av{iType},1),loC,hiC,'jet'),'filled');            
        scatter(x,y,30,mean(av{iType},1),'filled');
        colormap(myColorMap);
        colorbar('EastOutside');
        % set(gca,'Color','k');     % Use only if you don't plan to export these figures to illustrator!
        if(doCompareCells)          % If selectivity was measured
            indSelective = find(selectiveCells);
            for(iCell=1:length(indSelective))
                circle([x(indSelective(iCell)) y(indSelective(iCell))],3,10,'r-');
            end
        end
        drawnow;
        hold off;
    end
    supertitle('Amplitudes');
    % figure;
    % p = paired_comparisons([mean(av{1},1)' mean(av{2},1)'],{'f','c'}); title('Amplitudes');
    % text(1.5,max(mean(av{1},1)),sprintf('P = %s',myst(p(2,1))));    % Show p-value on the figure
end


%%% ------------ Processing of latencies (centers of mass over time)
centerlat = nan(nType,nCells);          % Latency centers of mass
for(iType = 1:nType)    
    temp = (1:size(av{iType},1))*av{iType}./sum(av{iType},1)*dt*1000;   % Center of mass, in ms
    temp = temp-mean(temp(~isnan(temp)));                               % Now unbias
    % temp = temp/std(temp(~isnan(temp)));                              % And maybe normalize
    temp(isnan(temp)) = 0;                                              % if the cell was silent, set centerLat to zero
    centerlat(iType,:) = temp;    
    
    temp = cumsum(av{iType},1);                                         % Now another approach, based on cumulative sum.
    for(iCell=1:nCells)
        lat = find(temp(:,iCell)>temp(end,iCell)*0.1,1,'first')*dt*1000;    % Also to ms
        if(isempty(lat))
            lat = NaN;
        end
        cumlat10(iType,iCell) = lat;
    end
end
if(doCompareCells)                                          % Everything that is about selectivity
    [rho,pval] = corr(kSelectiveFC',centerlat(cType,:)');   % Data needs to be columns for CORR to work
    p = polyfit(kSelectiveFC',centerlat(cType,:)',1);
    figure; hold on;
    plot(kSelectiveFC,centerlat(cType,:),'bo');    
    plot(get(gca,'XLim'),polyval(p,get(gca,'XLim')),'g-');
    xlabel('Selectivity for C against F'); ylabel('Weighted latency during responses to C');
    title(sprintf('R = %s, p = %s',myst(rho),myst(pval)));
    hold off;
    
    g = ~isnan(cumlat10(cType,:));
    [rho,pval] = corr(kSelectiveFC(g)',cumlat10(cType,g)');   % Data needs to be columns for CORR to work
    p = polyfit(kSelectiveFC(g)',cumlat10(cType,g)',1);
    figure; hold on;
    plot(kSelectiveFC,cumlat10(cType,:),'bo');    
    plot(get(gca,'XLim'),polyval(p,get(gca,'XLim')),'g-');
    xlabel('Selectivity for C against F'); ylabel('10% latency during responses to C');
    title(sprintf('R = %s, p = %s',myst(rho),myst(pval)));
    hold off;
    ylim([0 inf]);
end

% position_finder(x,y,cumlat10(cType,:))
% position_finder(x,y,centerlat(cType,:))

if(doLatencyMap)
    %loC = min(mean(av{1},1));          % For color scale
    %hiC = max(mean(av{1},1));
    figure;
    for(iType = 1:nType)
        if(nType>1)
            subplot(2,2,iType);
            title(key(iType));
            hold on;
        end        
        % scatter(x,y,15,colorscale(mean(av{iType},1),loC,hiC,'jet'),'filled');                    
        scatter(x,y,30,centerlat(iType,:),'filled');
        %colormap(myColorMap);
        colorbar('EastOutside');
        % set(gca,'Color','k'); % Use only if you don't plan to export these figures to illustrator!
        drawnow;
        hold off;
    end
    supertitle('Latencies');
    % figure;
    % p = paired_comparisons([mean(av{1},1)' mean(av{2},1)'],{'f','c'}); title('Amplitudes');
    % text(1.5,max(mean(av{1},1)),sprintf('P = %s',myst(p(2,1))));    % Show p-value on the figure    
end


% %%% ------------ Reference patch-clamped cell
% refCell = 47;   % 47, 46 and 37 are the only three cells in the vicinity (unless it was an invisible cell). Of these 47 is the most likely.
% dataE = abfload('C:\_Data\___Ca imaging\131213\2013_12_13_0001.abf');
% figure('Color','w');
% if(0)   % Two subplots nearby
%     subplot(1,2,1); hold on; title('Ephys');
%     for(q=1:12)
%         temp = squeeze(dataE(1:30000,1,q));
%         plot(temp-temp(1)+q);
%     end
%     hold off;
%     subplot(1,2,2); hold on; title('Ca');
%     for(q=1:12)
%         plot(S(q).dataS(:,refCell)*5+q);
%     end
%     hold off;
% else    % One plot with both types of data
%     eTime = (1:30000)/10;   % Ephys time, first 3 s in ms    
%     hold on; title(['Cell #' num2str(refCell)]);
%     timeShift = -520;    % in ms, ephys to the left
%     for(q=1:12)
%         temp = squeeze(dataE(1:30000,1,q));
%         plot(eTime+timeShift,temp-temp(1)+q,'b');             % Physiology   
%         n = size(S(q).dataS,1);     % n is different in different sweeps =/
%         cTime = (1:n)'/n*3000;      % Ca imaging time, same first 3 s.
%         plot(cTime,S(q).dataS(:,refCell)*5+q,'g');
%         xlim([0 2000]);
%     end
%     hold off;
%     drawnow;
% end



% %%% ------------ Amplitude variability
% timeSlice = 1:find(time>1000,1,'first');    % First second
% loC = min(std(av{1},1)./mean(av{1},1));     % For color scale
% hiC = max(std(av{1},1)./mean(av{1},1));
% figure;
% for(iType = 1:nType)
%     subplot(2,2,iType);
%     %scatter(x,y,15,colorscale(std(av{iType},1)./mean(av{iType},1),loC,hiC,'jet'),'filled','Marker','o');    
%     scatter(x,y,15,std(av{iType},1)./mean(av{iType},1),'filled','Marker','o');    
%     title(key(iType));
%     set(gca,'Color','k');  % Background color - use only if you don't plan to export these figures to illustrator!
%     drawnow;
% end
% supertitle('Response amplitude relative variation');
% % print(gcf,'-dmeta','-r600','C:\Users\Arseny\Documents\7_Ca imaging\figExport.emf') % EMF: does export, but AI doesn't read it
% % print(gcf,'-depsc2','-r600','C:\Users\Arseny\Documents\7_Ca imaging\figExport.eps') % EPS2: formally exports something, but it is cut in disgusting stripes
% % %%% NOTE: For some reason in this plot standard vector copying didn't work when background color was set to black. It worked without black background though.
% % %%%         At the same time for similar (almost identical?) figures below (latency & length) vector copying works even with black background. A mystery!
% 
% figure;
% p = paired_comparisons([var(av{1},1)' var(av{2},1)'],{'f','c'}); title('Relative variation');


% %%% Pseudo-Movie
% timeStep = 2;
% timeShift = [0 0 0 20];
% for(iType=[2 4]) % 2 = crash, 4 = back
%     figure('Color','w');
%     for(t=1:16)
%         subplot(4,4,t);        
%         scatter(x,y,7,colorscale(mean(av{iType}(timeShift(iType)+(t-1)*timeStep+(1:timeStep),:),1),0,0.2,'hot'),'filled');         
%         title(num2str(timeShift(iType)+(t-1)*timeStep));        
%         set(gca,'Color','k','XTick',[],'YTick',[]);
%         drawnow;
%     end
% end

% %%% -- Coordinates of enumerated cells, with a circle deliniating zones for comparison
% figure('Color','w');
% enum(x,y,1:length(x));
% hold on;
% circle([xc yc],circR);
% % plot(x(indIn),y(indIn),'.');    % Checking if it works
% hold off;
% set(gca,'DataAspectRatio',[1 1 1]);

% %%% --- In vs Out: Inside the circle vs outside
% figure('Color','w');
% for(iType=1:4);
%    subplot(4,1,iType);
%    hold on;
%    plot(time(time<1000),mean(av{iType}((time<1000), indIn),2),'r-');
%    plot(time(time<1000),mean(av{iType}((time<1000),~indIn),2),'b-');
%    legend({'In','Out'});
%    ylabel(stimName{iType},'FontSize',8);
%    set(gca,'FontSize',8);
%    hold off;
%    drawnow;
% end


% %%% --- Approximate response latencies, and lengthes (OLDER VERSION, PARTIALLY SUPERCEDED, NEEDS REWRITING)
% timeSlice = 1:find(time>1000,1,'first');    % First second
% loC = max(time(timeSlice(end)));  hiC = 0;  % For latencies color-coding. low_margin was assigned a high value, and high - a low value, so that they could reverse nicely.
% loClen = max(time(:)); hiClen = 1;          % For response length color-coding (also reverse-assignment).
% lat = [];   % Here latencies will be stored
% len = [];   % And here lengthes will be stored
% for(iType=1:nType)
%     for(iCell=1:m)        
%         % temp = find(av{iType}(timeSlice,iCell)==max(av{iType}(timeSlice,iCell)),1,'first'); % Find the maximum        
%         temp1 = find(cumsum(av{iType}(timeSlice,iCell),1)/sum(av{iType}(timeSlice,iCell))>0.1,1,'first'); % Based on cumulative sum threshold
%         temp2 = find(cumsum(av{iType}(timeSlice,iCell),1)/sum(av{iType}(timeSlice,iCell))>0.9,1,'first'); % Based on cumulative sum threshold        
%         lat{iType}(iCell) = time(temp1);                % Response latency
%         len{iType}(iCell) = time(temp2) - time(temp1);  % Response length
%     end
% %     lat{iType}(:) = lat{iType}(:) - min(lat{iType}(:));
%     loC = min(loC, min(lat{iType}(:)));
%     hiC = max(hiC, max(lat{iType}(:)));
%     loClen = min(loClen, min(len{iType}(:)));
%     hiClen = max(hiClen, max(len{iType}(:)));    
% end
% 
% figure; % -- Latencies
% colormap('jet');
% for(iType = 1:nType)
%     subplot(2,2,iType);    
%     scatter(x,y,10,lat{iType},'filled');
%     set(gca,'CLim',[loC hiC]);
%     title(key(iType));
%     set(gca,'Color','k');
%     colorbar('EastOutside','FontSize',8);    
%     drawnow;
% end
% supertitle('Response lantecy');
% 
% figure; % -- Lengthes
% colormap('jet');
% for(iType = 1:nType)
%     subplot(2,2,iType);    
%     scatter(x,y,10,len{iType},'filled');
%     set(gca,'CLim',[loClen hiClen]);
%     title(key(iType));
%     set(gca,'Color','k');
%     colorbar('EastOutside','FontSize',8);    
%     drawnow;
% end
% supertitle('Respones length');
% 
% figure; % -- Compare length and latency between F and C
% subplot(1,2,1);
% p = paired_comparisons([lat{1}(:) lat{2}(:)],{'f','c'}); title('Latencies');
% text(1.5,max(lat{1}),sprintf('P = %s',myst(p(2,1))));    % Show p-value on the figure
% subplot(1,2,2);
% p = paired_comparisons([len{1}(:) len{2}(:)],{'f','c'}); title('Lengthes');
% text(1.5,max(len{1}),sprintf('P = %s',myst(p(2,1))));    % Show p-value on the figure

% %%% ------------------ PCA --------------------------
% if(doPca)
%     timeSlice = 5:find(time>1000,1,'first');    % First second
%     amps = [];
%     trialsIncluded = 0;
%     for(iTrial = 1:nTrials)    
%         iType = mod(iTrial-1,nType)+1;    
%         if(iType==2)
%             amps = [amps S(iTrial).dataS(timeSlice,:)];        
%             trialsIncluded = trialsIncluded+1;
%         end
%     end
%     cellId = repmat(1:m,1,trialsIncluded);
%     [cPca,sPca,eigenvalues] = princomp(amps);
%     nComps = 4;
%     for(q=1:nComps)
%         c{q} = reshape(cPca(q,:),m,[]);
%     end
%     [~,i] = sort(mean(c{1},2));    % Cells sorted by amplitude
%     figure; myplot(c{1}(i,:));
%     ylabel('Cell number'); xlabel('Trial #');
%     title('Response intensity (1st component)');
% 
%     figure; plot(sPca(:,1:nComps));  % Component shapes
%     for(q=1:nComps); lab{q} = num2str(q); end;
%     legend(lab);
% 
%     map = [(31:-1:0)' zeros(32,1) (31:-1:0)'; zeros(32,1) (0:1:31)' zeros(32,1)]/31;
%     figure; colormap(map);
%     map = colormap;
%     for(q=1:nComps)
%         hA(q) = subplot(2,nComps,q);
%         scatter(x,y,10,mean(c{q},2)*eigenvalues(q),'filled');
%         %colorbar('EastOutside');
%         set(gca,'Color','k');
%         %set(gca,'Color',map(round(size(map,1)/2),:));      % Middle of the colormap
%         title(num2str(q));
%         subplot(2,nComps,nComps+q);
%         plot(sPca(:,q));
%     end
%     unifyc(hA);
% end

drawnow;



% %%% -------- Send data outside
% [b,a] = butter(3,1/2);
% useAverages = 0;
% sendSpikes = 0;
% if(useAverages) % Averages    
%     for(iType=1:nType)
%         res{iType} = filter(b,a,av{iType});
%     end
% else            % Full data (filtered)
%     for(iCell=1:nCells)
%         temp = [];
%         for(iTrial = 1:nTrials)
%             iType = mod(iTrial-1,nType)+1;
%             if(sendSpikes)
%                 temp = [temp S(iTrial).dataS(find(time>=latAdj(iType),1,'first')-1+ (1:ntimeticks),iCell)];
%             else                
%                 localData = S(iTrial).dataF(find(S(iTrial).time>=latAdj(iType)/1000,1,'first')-1+ (1:ntimeticks),iCell);                
%                 temp = [temp localData/localData(1)];  % Lots of backwards-compatibility (s, jerky timescale)
%             end
%         end
%         if(sendSpikes)
%             res{iCell} = max(0,filter(b,a,temp));
%         else
%             res{iCell} = temp;
%         end
%     end
% end

end


function unifyc(hA)
% Unify color scales
lo = []; hi = [];
for(q=1:length(hA))
    lim = get(hA(q),'CLim');
    if(isempty(lo)); lo = lim(1); else; lo = min(lo,lim(1)); end;
    if(isempty(hi)); hi = lim(2); else; hi = max(hi,lim(2)); end;
end
for(q=1:length(hA))
    set(hA(q),'CLim',[lo hi]);
end
end


function c = colorscale(x,lo,hi,name)
% I don't need this function anymore, but some functions above may still reference it.
% Better to get rid of it though. Jan 10 2014
x = x(:);
% x may be a horizontal vector, and I don't like it.
if(isempty(lo));     lo = min(x);   end
if(isempty(hi));     hi = max(x);   end
cmap = colormap(name);
n = size(cmap,1);
t = min(max(1,round(n*(x-lo)/(hi-lo))),n);
c = cmap(t,:);
end