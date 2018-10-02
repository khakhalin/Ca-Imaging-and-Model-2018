function shotgun_testing()
% shotgun_testing()
%
% A routine to test the effects of sparse sampling on network measurements

% Sep 15 2018: Created
% Oct 02 2018: Subsampling for large nCells removed. Some cleanup.


nSteps = 7;                         % How many rewire steps to perform
nCells = 500;                       % The modularity function works OK for up to 1000 cells, but then gets extremely slow at ~3000 (30 s for one calculation), 
                                    % and becomes completely unresponsive at 10000 (no result after several hours of waiting)
nRewire = floor(0.004*nCells^2);    % How many edges (not nodes!) to rewire at each run
nSample = 100;                      % How many cells to sample
nNetworks = 5;                      % How many different networks (independent experiments) to run
nSubsamples = 2;                    % How many substamples from each network to try
nCellsMax = 3000;                   % Max nCells possible (a safety thing, as when nCells is too high, everything just freezes)

useErrorBars = 0;
rewiringType = 'noise';             % options: maslov (degree preserving), random (move where the edges are pointing), noise (xor a bunch of edges)

hF = figure('Color','white');
listOfTitles = {'Modularity','Clustering','Hierarchy','Efficiency'};
listOfLimits = [1 0.4 0.2 1];   % X and Y upper limits for each of the measures

for(iPlot=1:4)                                                          % Create all axes with proper labels; we'll be platting into them
    ha(iPlot) = subplot(2,2,iPlot); 
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
    % A = create_modular_network(nCells,floor(nCells/10),0.9,0.1);	% Loosely connected rings - presumably cyclicity should be good, but not really
    A = create_modular_network(nCells,floor(nCells/10),0.5,0.2);	% A thin doughnut - decent modularity, low clustering
    
    for(iStep = 1:nSteps)
        fprintf('Network %d/%d - step %d/%d\n',iNetwork,nNetworks,iStep,nSteps);
        if(nCells>nCellsMax)
            error('Number of Cells is too high; in practice the program will never finish.');
        end
        ind = randperm(nCells);                                     % Shuffle indices
        ind = ind(1:nCells);                                        % Then take nSample first ones. That would be our sample
        B = A(ind,ind);                                             % Sub-network: Looking only at it, to make calculations feasible        
        
        [~,aMod] = modularity_dir(B);                               % Full network measuremens: modularity
        aClu = mean(clustering_coef_bd(B));                         % Clustering
        [aHie, ~] = network_flow(B);                                % Flow hierarchy    
        aEff = efficiency_wei(B);                                   % Global efficiency (average inverse shortest path between al pairs of nodes)

        iTry = 1;
        while(iTry<=nSubsamples)
            ind = randperm(nCells);                                 % Shuffle indices
            ind = ind(1:nSample);                                   % Then take nSample first ones. That would be our sample
            C = A(ind,ind);                                         % Sub-network
            try
                [~,yMod(iTry)] = modularity_dir(C);                 % Modularity sometimes produces an error (when the network gets disconnected?), so it's wrapped in a TRY
                yClu(iTry) = mean(clustering_coef_bd(C));           % Clustering
                [yHie(iTry), ~] = network_flow(C);                  % Flow hierarchy via Katz centrality                
                yEff(iTry) = efficiency_wei(C);                     % Efficiency
            catch                                                   % If it was a really bad try (for example, the graph become weirdly disconnected and something crashed)...
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
            plot(ha(3),aHie,yHie,'b.');
            plot(ha(4),aEff,yEff,'b.');
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