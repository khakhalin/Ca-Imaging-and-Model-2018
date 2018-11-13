function graph_structure_analyzer()
% graph_structure_analyzer()
%
% The function reads files with graph data one after another (the file names should be hard-coded in the header),
% analyzes their network properties, and then compares multiple measures across two groups. The 1/0 switches to turn
% different parts of analysis on and off are located below, as constants for the 1st function that the main function calls.
% All outputs go to figure and console; nothing is written in the external files.

% Aug 29 2017: Created
% Sep 01 2017: New network functions
% Sep 04 2017: Analysis of selectivity
% Sep 15 2017: Full week of improvement in fair edge thresholding between experiments
% Sep 29 2017: Edge analysis, spatial analysis
% Aug 09 2018: Revived.

onlyOneBrain = 0;       % When generating sample images, set this to 1, to process only a small subset of experiments (currently: 2 brains, see below).

useFullSet = 0;         % Set to 1 for full set (including noisy data); set to 0 for a clean subset (low noise, good graphs)
if(~onlyOneBrain)
    if(useFullSet)
        %%% --- stage 46 set
        folderName1 = 'C:\_Data\___Ca imaging\_caimg s46 results\';
        fileNames1 = {'140718a','140716b','140716a','140715','140711','140708b','140708a','140705a','140704b','140704a',...
            '140627','140626','140620','140619b','140619a','140613','140611','140610','140530b','140530a','140529','140528','140516','140502'};
        %%% --- stage 49 set
        folderName2 = 'C:\_Data\___Ca imaging\_caimg s49 results\';
        fileNames2 = {'140726','140724','140723','140722','140718b','140714','140710','140709','140707','140705b','140703',...
            '140612','140522','140521','140505','140408','140328','140326','140325','140318','140317','140314','140312','140311','140310'};
    else
        %%% --- stage 46 set
        folderName1 = 'C:\_Data\___Ca imaging\_caimg s46 results\';
        fileNames1 = {'140718a','140716b','140716a','140715','140711','140708b','140708a','140705a','140704b','140704a',...
            '140627','140626','140610','140516'};
        %%% --- stage 49 set
        folderName2 = 'C:\_Data\___Ca imaging\_caimg s49 results\';
        fileNames2 = {'140726','140724','140723','140722','140718b','140714','140710','140709','140707','140705b',...
            '140612','140521','140328','140314','140311','140310'};
    end
else %%% --- A minimal set for raw figures only
    folderName1 = 'C:\_Data\___Ca imaging\_caimg s46 results\';
    fileNames1 = {'140516'};
    %%% --- stage 49 set
    folderName2 = 'C:\_Data\___Ca imaging\_caimg s49 results\';
    fileNames2 = {'140310'};
end
    
n46 = length(fileNames1);   % How many s46s are there
n49 = length(fileNames2);   % How many s49s

data = [];                  % Place to collect network measurements, returned by the ANALYZE_GRAPH function below.
centralityBag = [];         % Place to collect info about all nodes
for(iExp=1:(n46+n49))
    if(iExp<=n46)
        fullName = [folderName1 'results-' fileNames1{iExp} '.mat'];
        stages(iExp) = 46;
    else
        fullName = [folderName2 'results-' fileNames2{iExp-n46} '.mat'];
        stages(iExp) = 49;
    end
    S = load(fullName);
    S = S.res;
    S.stage = stages(iExp);
    S.centralityBag = [];                           % To collect centrality data for every cell
    
    [newData,newLabels,centralityData] = analyze_graph(S,iExp);    % Make sure that both of these are rows, not columns
    
    centralityBag = [centralityBag; centralityData];
    labels = newLabels;
    data = [data; newData];
    
    %fprintf('\n');
end

if(~isempty(centralityBag))
    thisPathWithName = mfilename('fullpath');                               % Get the location of this file
    thisName = mfilename();
    localPath = thisPathWithName(1:(length(thisPathWithName)-length(thisName)));    % Path to the folder where this function lies
    
    fid = fopen([localPath 'sel_centrality_allcells.csv'],'w');              % csvwrite doesn't support headers, but we need a header
    fprintf(fid,'%s\n','nexp, sel, indegree, katz, spiking, insel');
    % [ones(nCells,1)*iExp, sel(:), inDegree(:), netrank, spiking(:), w'*sel(:)];
    fclose(fid);
    dlmwrite([localPath 'sel_centrality_allcells.csv'],centralityBag,'-append'); % write data to the end
end

fprintf('\n'); dispf(data);
if(~isempty(data))
    plotEverything(data,stages,labels);
end

end


function [netData,labels,centralityData] = analyze_graph(S,iExp)
sigLevel = 0.01;                % Significance level for accepting that the TE exists (doesn't matter if you use any of the adaptive techniques)
thresholdType = 'lax';          % Either 'strict' or 'lax' (recommended). With 'strict' many graphs end up empty or nearly empty.
wComboMode = 'mean';            % Three options of how wc, wf, and ws should be combined into one matrix: min, mean (recommended), max, c, f, s
minEdgeNumber = 0;              % How many edges in each graph to achieve. Set to 0 if it should not be adaptive. Recommended non-adaptive value: 20
maxEdgeNumber = 0;              % If set to above minEdgeNumber, shrinks large graphs down to this value by thresholding low ws out. Recommended: 50 (or better yet: 0)
minComponentSize = 0;           % Number of nodes in the largest weakly connected component. If set >1, supercedes minEdgeNumber.
forceAverageDegree = 1.0;       % if >0, instead of all other adaptive thresholding above, selects sigLevel that achieves this average degree. Recommended: 1.0
nShuffleNetworkMeasures = 0;    % How many times to reshuffle the matrix for all measurements. Set to 0 for no shuffling. Recommended: 50 or 100 (20 is too variable)
shuffleNetMode = 'degree';      % Two options: 'Erdos' (keeps n points and n edges, but reshuffle them), and 'degree' for degree-preserving (calls external function)

showEdgeNumbers = 0;            % Report cross-protocol replication stats, number of cells, edges, significance thresholds used etc.
showBasicStats = 0;             % Basic stats about response strength, correlations, etc.
showRawGraphs = 0;              % Set to 1 if graphs for each experiment need to be shown
doRawEdgeAnalysis = 0;          % Looking at raw edges, before thresholding
showEdgeAnalysis = 0;           % Analysis of weights, degrees, edge assymetry, and w<0 edges. The type of output has to be switched in the code block itself (sorry)
doPredictions = 0;              % Whether we want to look at stimulus identity prediction from brain activity
selToUse = 'FC';                % Which selectivity to use; options are: 'C' (combined), 'FC', and 'SC'
doSelectivity = 1;              % Where the selective cells are. Also calculates correlations with centrality measures and similar things
doSelectivityInGraph = 0;       % Analyze network properties and cell location. Relies on selectivity analysis

doConnectivityInSpace = 0;      % Spatial analysis for connectivity (for example, in which direction are edges facing, etc.)
nConnectivityShuffles = 100;    % To compare connectivity of original network to randomized network. Recommended: ~100
doNetworkMetrics = 0;           % Whether network properties need to be measured
doEnsembleConnectivity = 0;     % Load up ensembles identified by caimaging_pca.m , and compare connectivity within and between

auxFolder = 'C:\_Data\___Ca imaging\auxData\';  % Where the aux data is stored
D = load([auxFolder S.name '-auxData.mat']); D = D.D;

%%% Future outputs:
labels = {};                            % Labels for network measurements
netData = [];                           % Measurements themselves
centralityData = [];                    % Default, empty

% % --- Minimal cell distance: if it's too low for any of the trials, that's because there was ROI duplication in the dataset
% labels = [labels,'Min Cell Dist'];          
% minDist = max(S.xy(:,1))-min(S.xy(:,1));
% for(i=2:S.nCells)
%     for(j=1:i-1)
%         minDist = min(minDist, sqrt((S.xy(i,1)-S.xy(j,1))^2+(S.xy(i,2)-S.xy(j,2))^2));
%     end
% end
% data = [data minDist];

%%% --------------------------- Select which TE values will be analyzed
ncon = S.nCells*(S.nCells-1);               % Number of connections
S.teC.p = max(S.teC.p,eye(S.nCells));       % By mistake I set p-values for self-connections to a small value, not to 1. That's how it is now stored in results.
S.teF.p = max(S.teF.p,eye(S.nCells));       % Because of that, it needs to be explicitly fixed here. As self-connections don't exist, and so should have large p.
S.teS.p = max(S.teS.p,eye(S.nCells));


%%% ---------------- Raw edge weight anaysis
if(doRawEdgeAnalysis)
    rawTE = (S.teC.te + S.teF.te + S.teS.te)/3;
    % figure; hist(rawTE(:),100); drawnow();
    distSkew = skewness(rawTE(:),0);
    shareOnTop = sum(rawTE(:)>0.05)/length(rawTE(:));
    fprintf('%10s\t%8f\t%8f\n',S.name,distSkew,shareOnTop);
end


%%% --------- Find edges
if(minEdgeNumber>0 || minComponentSize>1 || forceAverageDegree>0)              
    sigLevel = 1e-12;                        % If we have to use adaptive threshold, start with the smallest p imaginable
end
flagNeedMoreEdges = 1;
while(flagNeedMoreEdges)
    switch wComboMode
        case 'c';   g = S.teC.p(:)<sigLevel;    % g stands for "good", as a variable to replace (:) with (g) in edge-wise calculations below
        case 'f';   g = S.teF.p(:)<sigLevel;
        case 's';   g = S.teS.p(:)<sigLevel;
        otherwise
            switch thresholdType
                case 'strict'
                    g = (S.teC.p(:)<sigLevel)&(S.teF.p(:)<sigLevel)&(S.teS.p(:)<sigLevel);  % Simple clear formula: all 3 should be significant, that's it
                case 'lax'
                    gp = S.teC.p(:).*S.teF.p(:).*S.teS.p(:);                                % Crazy formula that the product of 3 p-val should be < alpha^3
                    g = ( gp < sigLevel^3); 
            end
    end
    if(forceAverageDegree>0)                    % If we go with just setting a certain average degree (edge to node ratio)
        if(sum(g)/S.nCells>=forceAverageDegree)
          flagNeedMoreEdges = 0;
        else
            sigLevel = sigLevel*1.01;            % If not enough edges - increase p a little bit
            %gp = sort(gp); sigLevel = gp(find(gp>sigLevel,1,'first')); % Alternative approach that does not work
        end
    elseif(minComponentSize>1)                  % If minimal largest connected component size is set, use it for adaptive thresholding
        w = ones(size(S.teC.te)); 
        w(:) = w(:).*g(:);                                  % Prepare simplified adjacency matrix
        binComps = conncomp(digraph(w),'Type','weak');      % Find all weakly connected components of this graph
        largestComponent = sum(binComps==mode(binComps));   % Size of the largest connected component
        if(largestComponent>=minComponentSize)
            flagNeedMoreEdges = 0;
        else
            sigLevel = sigLevel*1.1;
        end
    else                                        % If minimal largest connected component size is not set, use number of edges
        if(sum(g)>=minEdgeNumber)
            flagNeedMoreEdges = 0;
        else
            sigLevel = sigLevel*1.1;
        end
    end
end

%%% --------- Using functional connectivity to infer structural connectivity
wc = S.teC.te; wc(:) = wc(:).*(g); wc(wc(:)<0) = 0; % Take tranfer entropies, remove those edges that aren't significant, and remove negative edges (if any).
wf = S.teF.te; wf(:) = wf(:).*(g); wf(wf(:)<0) = 0; % Note that w is a straight (unflipped) adjacency matrix, so G(w) makes sense.
ws = S.teS.te; ws(:) = ws(:).*(g); ws(ws(:)<0) = 0;

% Load all aux data. We need to load all in advance, as not all cells will go into the final analysis (see "core" below).
switch selToUse    
    case 'FC';  sel = D.selFC;              % Normal Cohen d for FC
    case 'SC';  sel = D.selSC;              % Same for SC
    case 'C';   sel = D.selC;               % Combo: (2C-S-F)/2 / (2SC+SS+SF)/4 
end
sel2 = D.selLogistic;                       % Different type of selectivigy (pseudo-r2 for logistic fit)
spiking = D.amp;                            % Total average spiking of each cell
corRaw = S.corRaw;                          % Raw correlation from the original data
corAdj = S.corAdj;                          % Adjusted correlation (effect of averages is partially taken out here)
ensembleInd = D.clustInd;                   % Load ensembles
originPoint = D.originPoint;                % Pre-calculated origin point for looming stimulus (a vector of [x y])

switch wComboMode
    case 'min';     w = min(wc,min(wf,ws));         % Leave only those connections that were obvious from all 3 stimuli
    case 'mean';    w = (wc+wf+ws)/3;               % Just a simple arithmetic mean; recommended.
    case 'max';     w = max(wc,max(wf,ws));
    case 'c';       w = wc;
    case 'f';       w = wf;
    case 's';       w = ws;
end

%%% --------- Shrinking the set based on W if requested
if(maxEdgeNumber>minEdgeNumber || forceAverageDegree>0)
    if(sum(w(:)>0)>maxEdgeNumber || sum(w(:)>0)/S.nCells>forceAverageDegree)
        threshold = min(w(:));
        maxWeight = max(w(:));
        largestComponent = maxEdgeNumber;                       % Just some arbitrary starting value
        while((sum(w(:)>0)>maxEdgeNumber && largestComponent>minComponentSize && maxEdgeNumber>0) || (sum(w(:)>0)/S.nCells>forceAverageDegree && forceAverageDegree>0))
            threshold = threshold + (maxWeight-threshold)*0.01; % Increase a bit
            wTemp = w;                                          % I want to try downsizing, but be able to back out of it if necessary
            wTemp(w(:)<threshold) = 0;
            if(minComponentSize==0)                             % if we are not playing with the largest component, just shrink it.
                w = wTemp;
            else                                                % If we are playing with largest component, which is more sensitive, be more careful:
                binComps = conncomp(digraph(wTemp),'Type','weak');  % Find all weakly connected components of this graph
                largestComponent = sum(binComps==mode(binComps));   % Size of the largest connected component            
                if(largestComponent>minComponentSize)               % Only go further if we haven't yet overshrunk w
                    w = wTemp;                                      % (and actually shrink it)
                end
            end
        end
    end
end

%%% --------- Reporting simple things
if(showEdgeNumbers)
    %%% Test: G = digraph([0 1 1 0 0 0; 0 1 0 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 0; 0 0 0 1 0 0; 0 0 0 0 0 0]); plot(G); conncomp(G,'Type','weak');    
    binComps = conncomp(digraph(w),'Type','weak');      % Identify and bin graph commponents
    nComponents = max(binComps);                        % How many components are there
    largestComponent = sum(binComps==mode(binComps));   % Largest component
    
    if(iExp==1)                                                 % If first brain, output column titles
        fprintf('Name        sigThresh   %%sigC   %%sigF   %%sigS   %%final   Ncon   Nfinal   H0     nCels   nPosE    nComp   LargestComp\n');
    end
    fprintf('%10s\t%5.2e\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',S.name,sigLevel,...
        sum(S.teC.p(:)<sigLevel)/ncon , sum(S.teF.p(:)<sigLevel)/ncon , sum(S.teS.p(:)<sigLevel)/ncon, sum(g)/ncon, ...
        ncon, sum(g) , round(sum(S.teC.p(:)<sigLevel)*sum(S.teF.p(:)<sigLevel)*sum(S.teS.p(:)<sigLevel)/ncon^2) , ...
        S.nCells, sum(w(:)>0), nComponents, largestComponent);
end

if(showBasicStats)
    if(sum(g)>0); corCF = corr(S.teC.te(g) , S.teF.te(g)); corCS = corr(S.teC.te(g) , S.teS.te(g)); 
    else;         corCF = NaN; corCS = NaN;    end
    mC = mean(S.teC.te(g)); sC = std(S.teC.te(g));
    mF = mean(S.teF.te(g)); sF = std(S.teC.te(g));
    mS = mean(S.teS.te(g)); sS = std(S.teC.te(g));
    fprintf('%10s\t%4.2f\t%4.2f\t%7.2e\t%7.2e\t%7.2e\t%7.2f\t%7.2f\n',S.name,...
        corCF , corCS, mC , mF , mS , (mC-mF)/sqrt((sC^2 + sF^2)/2) , (mC-mS)/sqrt((sC^2 + sS^2)/2));
end

%%% --------- Analysis of degree distributions
if(showEdgeAnalysis)
    G = digraph(w);
    indeg = histcounts(indegree(G),[0:10 20]-0.5);
    indeg = indeg/sum(indeg);
    oudeg = histcounts(outdegree(G),[0:10 20]-0.5);
    oudeg = oudeg/sum(oudeg);
    switch 'oudeg'                                      % Choose what to report
        case 'basic'
            if(iExp==1); fprintf('Name    %%symmetric     %%w<0\n'); end
            fprintf('%10s\t%10f\t%10f\n', S.name, sum(sym(:))/sum(w(:)), sum(inhib(:))/sum(w(:)));
        case 'indeg'
            %fprintf('%10s\t', S.name);
            if(iExp==1); fprintf('Reporting in-degrees:\n'); end
            dispf([indeg std(indegree(G))],'%5.2f');
        case 'oudeg'
            %fprintf('%10s\t', S.name);
            if(iExp==1); fprintf('Reporting out-degrees:\n'); end
            dispf([oudeg std(outdegree(G))],'%5.2f'); 
        case 'weights'
            if(iExp==1); fprintf('Showing (but not reporting) histograms of weights...\n'); end
            figure; histogram(w(:)); set(gca,'YScale','log'); xlim([0 0.1]); drawnow();
            %%% Pretty useless overall, so just look at these figures and realize it
    end
end

%%% --------- Get rid of unconnected nodes and tidy up the data
[w,indCore] = network_core(w);                         	% Leave only the core of the graph
w = w;                                                  % We are not flipping the matrix here, as all TEs are stored in direct, unflipped form. G(w) is a correct graph.
g = w(:)>0;                                             % Good edges for the purposes below
sel = sel(indCore==1);                                  % Selectivity
sel2 = sel2(indCore==1);                                % Logistic-fit-based selectivity
spiking = spiking(indCore==1);                                  % Oveerall amplitudes
corRaw = corRaw(indCore,indCore);                       % Raw correlations between spike-trains
corAdj = corAdj(indCore,indCore);                       % Correlations between spike-trains with average response for each cell taken out
nCells = size(w,1);                                     % New nCells as of now
amps = S.totalResponses(:,indCore==1);                  % Total response amplitudes for every neuron: nStim by nCells
ensembleInd = ensembleInd(indCore==1);

nStim = size(amps,1);
stimType = mod((1:nStim)-1,3)+1;                        % 1=c, 2=f, 3=s

%%% --------- Analysis of individual edges
if(showEdgeAnalysis)
    sym = w.*w';                                        % Symmetrical (bidirectional) connections    
    inhib = w.*(corAdj<0);                              % Not-excitatory connections (inhibitory or inactivating)
    %%% The reporting should be here, but it got killed by recent reorganization of the code. Add it back if necessary
end


if(doPredictions) %%% ----- How well can we predict stimulus identity from the response
    % Most useful variable: amps = total response, nStim by nCells    
    actuallyCrash = (stimType==1)';                                                  	% In experimental files, 1=C, as it gows CFS    
    firstHalf = [ones(floor(nStim/2),1); zeros(ceil(nStim/2),1)];                       % First half of all responses    
    warning('off');
    b = mnrfit(amps(firstHalf==1,:),actuallyCrash(firstHalf==1)+1,'model','nominal');   % Logistic regression on training half. +1 is needed to turn 0/1 into 1/2 (categories)
    warning('on');
    prediction = mnrval(b,amps(firstHalf~=1,:),'model','nominal');                      % Prediction on the testing half of the dataset
    predictionOfCrash = prediction(:,2)>prediction(:,1);
    predictionQuality = (sum(predictionOfCrash & actuallyCrash(firstHalf~=1))/sum(actuallyCrash(firstHalf~=1)) + ...
        sum(~predictionOfCrash & ~actuallyCrash(firstHalf~=1))/sum(~actuallyCrash(firstHalf~=1)))/2; % FPR+FNR/2 = expected accuracy on a balanced 50/50 set    
    if(iExp==1)                                                 % If first brain, output column titles
    	fprintf('  Name      prediction\n');
    end
    fprintf('%7s\t%8.2f\n',S.name,predictionQuality);
end

    
if(doSelectivity)
    if(0)   % Troubleshooting
        [froms(w(:)>0)'; tos(g)'; sel2(froms(g)); sel2(tos(g)); sel2(froms(g))<sel2(tos(g))]'
        figure('Color','white'); 
        subplot(1,2,1); plot(digraph(w),'NodeLabel',{},'NodeCData',sel ); set(gca,'visible','off'); drawnow();
        %subplot(1,2,2); plot(digraph(w),'NodeLabel',{},'NodeCData',sel2); set(gca,'visible','off'); drawnow();
        subplot(1,2,2); plot(digraph(w),'NodeCData',sel2); set(gca,'visible','off'); drawnow();
        error();
    end
    froms = repmat((1:nCells)',1,nCells);   % Different from-neurons are different rows, so we create a column of different rows and replicate right
    tos = froms';
    
    compRowLog = sel2(froms(g))<sel2(tos(g));
    pBiasLog = signrank(sel2(froms(g)),sel2(tos(g)));       % 
    shareSimpleLog = sum(compRowLog)/length(compRowLog);
    shareFancyLog  = compRowLog*w(g)/sum(w(g));
    
    G = digraph(w);
    if(0)   % Only analyze largest connected component. As it turns out, does not improve anything
        binComps = conncomp(G);
        nodeIds = find(binComps==mode(binComps));
        G = subgraph(G,nodeIds);
    else    % Analyze full graph
        nodeIds=(1:nCells);
    end
    
    inDegree = sum(w,1);
    ouDegree = sum(2,1);
    netrank = myCentrality(w,'netrank');                                    % Katz centrality
    clu = myCentrality(w,'clustering');                                     % Clustering coeff for every node
    %%% --- Correlations on nodes
    rSelDeg = corr(sel(:),inDegree(:));                                     % Correlation selectivity to in degree
    rSelIn2Ou = corr(sel(:),myCentrality(w,'gatherer'));                    % Selectivity to "almost ratio of inputs to outputs"
    rSelKatz  = corr(sel(:),netrank);                                       % Katz centrality: whether the cell is an integrator
    rSelRevKatz  = corr(sel(:),myCentrality(w,'revnetrank'));               % Reverse Katz centrality: whether the cells is an influencer
    rSelClu  = corr(sel(:),clu);                                            % Sel vs Clustering coefficient for this node
    rAmpClu = corr(spiking(:),clu);                                             % Response amp vs clustering coeff
    rAmpSel = corr(spiking(:),sel(:));                                          % Are selective cells more spiky?
    
    centralityData = [ones(nCells,1)*iExp, sel(:), inDegree(:), netrank, spiking(:), w'*sel(:)];
    % The last term here is the average selectivity of cells connected to this cell
        
    %%% --- Correlations on edges
    % selAssort = sel(:)'*w*sel(:)/sum(w(:))-((sum(sel*w)+sum(w*sel'))/(2*sum(w(:))))^2;  % Are selective cells are more likely to be connected? OBSOLETE, see below
    % Note that model_stdp_tester has the same formula but a different understanding of w (it is flipped there).
    % However the formula works either way, as (a'wa)' = a'w'a. And as it's a scalar, it means that w' or w doesn't matter here.    
    shi = repmat(sel(:)',nCells,1);                 % Matrix of properties for sending neurons
    sho = shi';                                     % Matrix of properties for receiving neurons
    selAssort = weightedcorr(shi(:),sho(:),w(:));
    
    selComparisonMatrix = bsxfun(@plus,sel(:)',-sel(:));                        % Matrix of sel growth: sel_j - sel_i (to interact with w, where w_ij is an edge from i to j)
    wTh = quantile(w(:),0.75);                                                  % Threshold for "strong edges"
    shESelGrow = sum((selComparisonMatrix(:)>0).*(w(:)>wTh))/sum(w(:)>wTh);     % What share of strong-ish edges increases Sel?
    selEGrowth = sum(selComparisonMatrix(:).*w(:))/sum(w(:));                   % Average weighted change in selectivity over an edge
    
    %%% --------- Degree distribution analysis -----
    temp = sort(w(:));                  % Sort all weights
    edgeThreshold = temp(end-nCells);   % find a thershold so that about nCells weights would pass it (similar to how we do it in experiments, to make E=V)
    [histDegIn,~] = histcounts(sum(w>edgeThreshold,2),(0:10)-0.5);	% In-degree distributions (of a thresholded graph)
    [histDegOu,~] = histcounts(sum(w>edgeThreshold,1),(0:10)-0.5);	% Out-degree distributions. Could have just used round() here.
    histDegIn = max(histDegIn,1)/nCells;                % Limit above zero to make log() possible
    histDegOu = max(histDegOu,1)/nCells;
    x = log(1:8);                                       % For fitting, skip deg==0 (arguably deg=1 could have been skipped as well, but I decided to include it)
    degFitIn = polyfit(x,log(histDegIn(2:9)),1);        % Fit of in-degrees, from 1 to 8
    gammaIn = degFitIn(1);
    degFitOu = polyfit(x,log(histDegOu(2:9)),1);        % Fit of out-degrees
    gammaOu = degFitOu(1);     
    
    %figure(hF); plot(sel(:),sel2(:),'.'); xlabel('Selectivity for C'); ylabel('abs(pseudo-R2) for log-regression of C'); drawnow();
    %figure(hF); plot(flowNode(:),sel2(nodeIds)','.'); xlabel('Flow to node'); ylabel('abs(pseudo-R2) for log-regression of C'); drawnow();
    if(1) % Troubleshooting figure
        hF = findobj('type','figure','name','Selectivity');
        if(isempty(hF)); hF = figure('Color','white','name','Selectivity'); hold on; end;
        set(gca,'XScale','log','Yscale','log');
        %plot(netrank,sel(nodeIds)','.'); 
        plot(0:9,histDegIn,'b.-');
        plot(0:9,histDegOu,'ro-');
        %xlabel('Cute measure'); ylabel('Selectivity'); 
        drawnow();    
    end

    if(iExp==1)                                                 % If first brain, output column titles
    	fprintf('  Name      rSelInDeg  rSelIntoOut rSelKatz    rSelRevKatz  rSelClu     rAmpClu     selAssort   shESelGrow  selEGrowth  gamma     rAmpSel\n');
    end
    fprintf('%7s\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\n',S.name,...
        rSelDeg,rSelIn2Ou,rSelKatz,rSelRevKatz,rSelClu,rAmpClu,selAssort,shESelGrow,selEGrowth,(gammaIn+gammaOu)/2,rAmpSel);
end


if(showRawGraphs)
    if(exist('sel2','var'))
        figure('Color','white','Name',S.name); plot(digraph(w),'NodeLabel',{},'NodeCData',sel2); set(gca,'visible','off'); drawnow();
        figure('Color','white','Name',[S.name 'raw']); plot(digraph(w),'NodeLabel',{},'NodeCData',sel2,'XData',S.xy(indCore,1),'YData',S.xy(indCore,2)); set(gca,'visible','off'); drawnow();
    else
        figure('Color','white','Name',S.name); plot(digraph(w),'NodeLabel',{},'NodeCData',sel); set(gca,'visible','off'); drawnow();
    end
end


%%% --------- Analysis of spatial arrangement of selectivity  and connectivity
x = S.xy(indCore,1);
y = S.xy(indCore,2);
if(doSelectivityInGraph)    
    distOriginToNode = sqrt((x(:)-originPoint(1)).^2 + (y(:)-originPoint(2)).^2);
    [rSelDist,pval] = corr(distOriginToNode(:),sel(:));                     % Is distance from the center correlated with selectivity?

    if(iExp==1); fprintf('Name       rSelDist       p_SelDist\n'); end
    fprintf('%8s\t%8.2f\t%8s\n',S.name,rSelDist,myst(pval));
    
    hF = findobj('type','figure','name','Spatial');
    if(isempty(hF)); 
        hF = figure('Color','white','name','Spatial'); 
        h.a1 = subplot(2,2,1); hold on; colormap(brewermap(20,'OrRd')); title('s45-46'); xlim([0 120]); ylim([0 120]); set(h.a1,'Ydir','reverse');
        h.a2 = subplot(2,2,2); hold on; colormap(brewermap(20,'OrRd')); title('s48-49'); xlim([0 120]); ylim([0 120]); set(h.a2,'Ydir','reverse');
        h.a3 = subplot(2,2,3); hold on; xlabel('Distance'); ylabel('Selectivity'); xlim([0 120]); %ylim([0 1]);
        h.a4 = subplot(2,2,4); hold on; xlabel('Distance'); ylabel('Selectivity'); xlim([0 120]); %ylim([0 1]);
        set(hF,'UserData',h);
    else
        h = get(hF,'UserData');
    end
    if(S.stage<47)  % Young - left plot
        scatter(h.a1,x,y,40,sel','filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','none');    
        plot(h.a3,distOriginToNode,sel,'b.');
    else
        scatter(h.a2,x,y,40,sel','filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','none');    
        plot(h.a4,distOriginToNode,sel,'b.');
    end
    drawnow();
end

if(doConnectivityInSpace)    
    if(iExp==1); fprintf('Calculation started...\n'); end;    
    xmat = repmat(x(:),1,nCells);                               % Repeat each position nCell times ot calculate pairwise differences
    ymat = repmat(y(:),1,nCells);
    pwd = sqrt((xmat-xmat').^2 + (ymat-ymat').^2);              % Matrix of pairwise distances
    avDistConn = sum(w(:).*pwd(:))/sum(w(:));                   % Average distance on current matrix
    avDistRand = sum(pwd(:))/(nCells*(nCells-1));               % Avereage distances if all were cocnnected to all, excluding self-connections
    xd = xmat'-xmat; yd = ymat'-ymat;                           % Vector from every cell to every other cell: [xd yd]. xd_ij is from i to j
    
    %%% At this point we should have originPoint var that was read from the D structure (well above).
    vectFromOriginX = (xmat'+xmat)/2-originPoint(1);            % Vector to the center of the edge connecting two nodes
    vectToOriginY = (ymat'+ymat)/2-originPoint(2);
    distFromOrigin = sqrt(vectFromOriginX.^2 + vectToOriginY.^2);
    origDotProd = (xd(:).*vectFromOriginX(:) + yd(:).*vectToOriginY(:)).*w(:)./pwd(:)./distFromOrigin(:)/sum(w(:));	% Normalized dot-products: cos(alpha)
    shareAway = sum(origDotProd>0)/sum(w(:)>0);                 % Share of edges that face away rather than towards the origin
    fromOrigin = sum(origDotProd(~isnan(origDotProd)));         % As we don't filter w_ii there are some 0 divided by 0s there. Dirty coding, sorry :\
    
    distOriginToNode = sqrt((xmat-originPoint(1)).^2 + (ymat-originPoint(2)).^2);
    distChange = distOriginToNode' - distOriginToNode;
    avChangeInDistance = sum(distChange(:).*w(:))/sum(w(:));    
    
    if(0) % Visual test that co-alignment is calculated correctly
        figure; hold on;
        plot([xmat(g) xmat(g)+xd(g)]',[ymat(g) ymat(g)+yd(g)]','-','Color',[1 1 1]*0.9); % Sketch edges
        plot(x,y,'k.');                             % nodes
        plot(originPoint(1),originPoint(2),'ks');   % sketched edges
        for(i=1:nCells)
            for(j=1:nCells)
                k = (j-1)*nCells + i;               % place within the xmat-like matrix
                if(w(i,j)>0)
                    if(distChange(k)>0)             % alternatively put origDotProd(k) here
                        plot([x(i) x(i)+(x(j)-x(i))*0.9],[y(i) y(i)+(y(j)-y(i))*0.9],'r-');
                    else
                        plot([x(i) x(i)+(x(j)-x(i))*0.9],[y(i) y(i)+(y(j)-y(i))*0.9],'b-');
                    end
                end
            end
        end
        drawnow; 
    end    
    
    avDistShfl = [];
    fromCenterShfl = [];
    shareAwayShfl = [];
    wShfl = zeros(size(w));
    for(iTrial=1:nConnectivityShuffles)                         % The same measure on a randomly rewired network (degree-preserving)
        wShfl = network_rewire(w);                              % Maslov rewiring                
        avDistShfl(iTrial) = sum(wShfl(:).*pwd(:))/sum(wShfl(:)); % Reshuffles connections, but keeps biases due to hub cells         
        temp = (xd(:).*vectFromOriginX(:) + yd(:).*vectToOriginY(:)).*wShfl(:)./pwd(:)./distFromOrigin(:)/sum(w(:));	% Normalized dot-product: cos(alpha)
        shareAwayShfl(iTrial) = sum(temp>0)/sum(wShfl(:)>0);
        fromOriginShfl(iTrial) = sum(temp(~isnan(temp)));
        avChangeInDistanceShfl(iTrial) = sum(distChange(:).*wShfl(:))/sum(wShfl(:));
    end
    [~,pConnectedAreCloser] = ttest(avDistShfl-avDistConn);    
    [~,pEdgesPointAway] = ttest(fromOriginShfl-fromOrigin);
    [~,pChangeInDistance] = ttest(avChangeInDistanceShfl-avChangeInDistance);
    
    if(iExp==1)                                                 % If first brain, output column titles        
        fprintf('  Name         avDistConn  avDistSfhl  avDistRand  P-dist      avAway      avAwayShufl P-Away       shareAway   shAwShfl   distCh      dChShfl    P-distCh\n');
    end    
    fprintf('%8s\t%8.2f\t%8.2f\t%8.2f\t%8s\t%8.2f\t%8.2f\t%8s\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8s\n',S.name,avDistConn,mean(avDistShfl),avDistRand,myst(pConnectedAreCloser),...
        fromOrigin,mean(fromOriginShfl),myst(pEdgesPointAway),shareAway,mean(shareAwayShfl),avChangeInDistance,mean(avChangeInDistanceShfl),myst(pChangeInDistance));
end


%%% --------- Measure network properties, with reshuffling (if needed)
if(doNetworkMetrics)                         
    if(nShuffleNetworkMeasures==0); nShuffleNetworkMeasures=1; end           % This would also be no reshuffling    
    fprintf('%8s',S.name);                      % No reporting here, as this data is returned to the higher level function, where it is processed
    for(is=1:nShuffleNetworkMeasures)           % if shuffling is not needed, just set nShuffle at the very top to 0        
        if(nShuffleNetworkMeasures>1)           % If more than one run, we reshuffle
            fprintf('.');                       % Cheap human's progress line
            switch(shuffleNetMode)
                case 'Erdos'
                    v = w;
                    v(:) = v(randperm(length(v(:))));       % Just completely reshuffle the edges
                case 'degree'
                    v = network_rewire(w,length(w(:))*3);	% Degree-preserving rewiring; rewires about 98% of all pairwise connections. That's a slow part.
                otherwise
                    error('Wrong type of shuffling: correct shuffleNetMode variable');
            end
        else
            v = w;                              % If only 1 run, just do the actual thing
        end
        newLine = [];
        
        temp = efficiency_wei(v);            newLine = [newLine temp];          labels = [labels, 'Efficiency'];  % --- Global efficiency        
        temp = assortativity_wei(v,1);       newLine = [newLine temp];          labels = [labels, 'Ou-In Assrt']; % 1 out-in, 2 in-out, 3 out-out, 4 in-in
        temp = assortativity_wei(v,2);       newLine = [newLine temp];          labels = [labels, 'In-Ou Assrt'];
        temp = assortativity_wei(v,3);       newLine = [newLine temp];          labels = [labels, 'Ou-Ou Assrt'];
        temp = assortativity_wei(v,4);       newLine = [newLine temp];          labels = [labels, 'In-In Assrt'];
        [~,temp] = modularity_dir(v);        newLine = [newLine temp];          labels = [labels, 'Modularity'];
        temp = clustering_coef_wd(v); temp = mean(temp(~isinf(temp)));  newLine = [newLine temp];          labels = [labels, 'Clustering'];        
        [flow, revFlow] = network_flow(v);   newLine = [newLine flow, revFlow]; labels = [labels, 'Flow', 'RevFlow'];
        temp = myCyclicity(v);               newLine = [newLine temp];          labels = [labels, 'Cyclicity'];
        
        netData = [netData; newLine];
    end
    fprintf('\n');                          % New row (we are printing names out, and maybe some dots, just to entertain the user)
    
    if(nShuffleNetworkMeasures>1)
        netData = mean(netData(2:end,:));   % Calculate the average
    end
end

%%% ----------------------------- Connectivity within ensembles --------------------
if(doEnsembleConnectivity)
    nClusters = max(ensembleInd);
    w_within = []; w_between = [];          % Weight of edges connecting within and between clusters    
    for(iCluster=2:nClusters)               % Check whether ensembles are closer to each other than between then, in space, and on a graph
        temp = w(ensembleInd==iCluster,ensembleInd==iCluster);
        w_within = [w_within; temp(:)];
        for(jCluster=1:(iCluster-1))
            temp = w(ensembleInd==iCluster,ensembleInd==jCluster);    % w is not flipped here, but it doesn't matter as we combine both directions
            w_between = [w_between; temp(:)];
            temp = w(ensembleInd==jCluster,ensembleInd==iCluster);
            w_between = [w_between; temp(:)];
        end
    end    
    [~,pval] = ttest2(w_within,w_between);
    % fprintf('N clusters: %2d; Within: %5.3f; between: %5.3f (p=%s); var explained: %4.2f\n',...
    %     nClusters,mean(bagWithin),mean(bagBetween),myst(pval),varExplained); % Console reporting
    withinClusterPreference = mean(w_within)/mean(w_between);
    
    if(iExp==1)                                                 % If first brain, output column titles       
        fprintf('  Name \n');
    end    
    fprintf('%8s\t%8.2f\t%8s\n',S.name,withinClusterPreference,myst(pval));
end


end



function plotEverything(data,key,labels)
% Prepares data for every column and calls "compare_columns" on it.

key = key(:);
nVar = size(data,2);
figure('Color','white');
plotn = 2; %ceil(sqrt(nVar)/1.5);
plotm = ceil(nVar/2); %ceil(nVar/plotn);
for(q=1:nVar)
    subplot(plotn,plotm,q);
    k = unique(key);
    t = [];    
    for(ik=1:length(k))
        t = concatenan(t,data(key==k(ik),q));
    end
    pval = compare_columns(t,k);
    title(sprintf('p=%s',myst(pval(1,2))),'FontWeight','Normal');
    ylabel(labels{q});
    set(gca,'FontSize',8);
end

end

