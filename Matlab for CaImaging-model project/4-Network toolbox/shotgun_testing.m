function shotgun_testing()
% shotgun_testing()
%
% A routine to test the effects of sparse sampling on network measurements

nSteps = 7;                         % How many rewire steps to perform
nCells = 500;                      % The modularity function works OK for up to 1000 cells, but then gets extremely slow at ~3000 (30 s for one calculation), 
                                    % and becomes completely unresponsive at 10000 (no result after several hours of waiting)
nCellsMax = 1000;                   % If nCells is greater than nCellsMax, we start subsampling estimations for some functions (modularity)
nRewire = floor(0.004*nCells^2);    % How many edges (not nodes!) to rewire at each run
nSample = 100;                      % How many cells to sample
nNetworks = 5;                      % How many different networks (independent experiments) to run
nSubsamples = 2;                    % How many substamples from each network to try

useErrorBars = 0;
rewiringType = 'noise';             % options: maslov (degree preserving), random (move where the edges are pointing), noise (xor a bunch of edges)

hF = figure('Color','white');
listOfTitles = {'Modularity','Clustering','Reach','Rev Reach','Efficiency'};
listOfLimits = [1 0.4 0.2 0.2 1];   % X and Y upper limits for each of the measures
for(iPlot=1:6)
    ha(iPlot) = subplot(2,3,iPlot); 
    if(iPlot<=length(listOfTitles))
        title(listOfTitles{iPlot});
    end
    hold on;
    if(iPlot<=length(listOfLimits))                                     % Cosmetics, to make scale at different plots different
        limit = listOfLimits(iPlot);
    else
        limit = 1;
    end
    plot([0 limit],[0 limit],':','Color',[1 1 1]*0.8);                  % Diagonal
    xlim([0 limit]); ylim([0 limit]);
end
   
for(iNetwork = 1:nNetworks)    
    %%% Create a fake network

    % A = create_modular_network(nCells,floor(nCells/10),0.7,0.1);	% DA string of connected clumps - good modularity, high clustering
    A = create_modular_network(nCells,floor(nCells/10),0.5,0.2);	% A thin doughnut - decent modularity, low clustering
    % A = create_modular_network(nCells,floor(nCells/10),0.9,0.1);	% Loosely connected rings - presumably cyclicity should be good, but not really
    for(iStep = 1:nSteps)
        fprintf('Network %d/%d - step %d/%d\n',iNetwork,nNetworks,iStep,nSteps);
        ind = randperm(nCells);                                     % Shuffle indices
        ind = ind(1:min(nCellsMax,nCells));                         % Then take nSample first ones. That would be our sample
        B = A(ind,ind);                                             % Sub-network: Looking only at it, to make calculations feasible        
        
        [~,aMod] = modularity_dir(B);                               % Full network measuremens: modularity
        aClu = mean(clustering_coef_bd(B));                         % Clustering
        [aFlo, aRfl] = network_flow(B);                             % Reach via Katz centrality        
        aEff = efficiency_wei(B);                                   % Global efficiency (average inverse shortest path between al pairs of nodes)

        iTry = 1;
        while(iTry<=nSubsamples)
            ind = randperm(nCells);                                 % Shuffle indices
            ind = ind(1:nSample);                                   % Then take nSample first ones. That would be our sample
            C = A(ind,ind);                                         % Sub-network
            try
                [~,yMod(iTry)] = modularity_dir(C);                 % Modularity sometimes produces an error (when the network gets disconnected?), so it's wrapped in a TRY
                yClu(iTry) = mean(clustering_coef_bd(C));           % Clustering
                [yFlo(iTry), yRfl(iTry)] = network_flow(C);         % Reach via Katz centrality
                temp = 1./C; temp(isnan(temp)) = 0;
                yEff(iTry) = efficiency_wei(temp);                  % Efficiency
            catch                                                   % If it was a bad try (for example, the graph become weirdly disconnected and something crashed)...
                continue                                            % ...just try again
            end
            iTry = iTry+1;                                          % Otherwise, if the results were fine, remember them
        end
        
        if(useErrorBars)
            bMod = mean(yMod);
            sMod = std(yMod);
            errorbar(ha(1),aMod,bMod,sMod,'b-','CapSize',0);     plot(ha(1),aMod,bMod,'b.');
        else
            plot(ha(1),aMod,yMod,'b.');
            plot(ha(2),aClu,yClu,'b.');
            plot(ha(3),aFlo,yFlo,'b.');
            plot(ha(4),aRfl,yRfl,'b.');
            plot(ha(5),aEff,yEff,'b.');
        end
        drawnow();

        if(iStep<nSteps)                                                % Prepare for next run: time to rewire the big matrix
            switch(rewiringType)
                case 'maslov'
                    A = network_rewire(A,nRewire);                      % Degree-preserving (too conservative)
                case 'random'
                    ind = randperm(nCells);                             % Wild rewiring
                    ind1 = ind(1:nRewire);
                    ind2 = [ind(2:nRewire) ind(1)];
                    A(ind1,ind1) = A(ind2,ind2);
                case 'noise'
                    A = xor(A,floor(rand(nCells)+(nRewire/nCells^2))^iStep); % Adaptive rewiring speed, to make the values a bit more uniformly distributed
            end
        end
    end
end

end