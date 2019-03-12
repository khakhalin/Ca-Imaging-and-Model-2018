function caimaging_structure(whichSet)
% caimaging_structure()
%
% Analysis of network structure.

% Mar 20 2017: Created.
% Aug 25 2017: Updated for TE calculation, fixed issues with shorter trials.
% Sep 04 2017: Now it should update data without overwriting TEs
% Sep 16 2017: Reworked some of the low-level analysis (doesn't affect TE calculations)
% Oct 05 2017: As of today PCA and mugshot functionality is considered abandoned here, and migrated to caimaging_pca

res = [];
iFolder = 0;

if(nargin<1); whichSet = 1; end;

if(whichSet==1)
    %%% --- stage 49 set
    folderBaseIn = 'C:\_Data\___Ca imaging\_caimg s49 mat\';
    folderBaseOut = 'C:\_Data\___Ca imaging\_caimg s49 results\';

    iFolder = iFolder+1; folderName{iFolder} = '140726'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140724'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140723'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140722'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140718b'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140714'; % done. Late responses (rebounds) for flash.
    iFolder = iFolder+1; folderName{iFolder} = '140710'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140709'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140707'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140705b'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140703'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140612'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140522'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140521'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140505'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140408'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140328'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140326'; % done. cfs2 is shorter than others
    iFolder = iFolder+1; folderName{iFolder} = '140325'; % done. Noisy
    iFolder = iFolder+1; folderName{iFolder} = '140318'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140317'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140314'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140312'; % done
    iFolder = iFolder+1; folderName{iFolder} = '140311'; % done. Noisy, except for 2 cells, but produced a decent graph
    iFolder = iFolder+1; folderName{iFolder} = '140310'; % done. Short recording
else
    %%% --- stage 46 set
    folderBaseIn = 'C:\_Data\___Ca imaging\_caimg s46 mat\';
    folderBaseOut = 'C:\_Data\___Ca imaging\_caimg s46 results\';
    
    iFolder = iFolder+1; folderName{iFolder} = '140718a';
    iFolder = iFolder+1; folderName{iFolder} = '140716b';
    iFolder = iFolder+1; folderName{iFolder} = '140716a';
    iFolder = iFolder+1; folderName{iFolder} = '140715';
    iFolder = iFolder+1; folderName{iFolder} = '140711';
    iFolder = iFolder+1; folderName{iFolder} = '140708b';
    iFolder = iFolder+1; folderName{iFolder} = '140708a';
    iFolder = iFolder+1; folderName{iFolder} = '140705a';
    iFolder = iFolder+1; folderName{iFolder} = '140704b';
    iFolder = iFolder+1; folderName{iFolder} = '140704a';
    iFolder = iFolder+1; folderName{iFolder} = '140627';
    iFolder = iFolder+1; folderName{iFolder} = '140626'; % Weak file, used for some testing
    iFolder = iFolder+1; folderName{iFolder} = '140620';
    iFolder = iFolder+1; folderName{iFolder} = '140619b';
    iFolder = iFolder+1; folderName{iFolder} = '140619a';
    iFolder = iFolder+1; folderName{iFolder} = '140613';
    iFolder = iFolder+1; folderName{iFolder} = '140611';
    iFolder = iFolder+1; folderName{iFolder} = '140610';
    iFolder = iFolder+1; folderName{iFolder} = '140530b';
    iFolder = iFolder+1; folderName{iFolder} = '140530a';
    iFolder = iFolder+1; folderName{iFolder} = '140529';
    iFolder = iFolder+1; folderName{iFolder} = '140528';
    iFolder = iFolder+1; folderName{iFolder} = '140516';
    iFolder = iFolder+1; folderName{iFolder} = '140502';
end

%%% ---------------- Control flags ----------------

flagOneBrainOnly = 1;                   % Set to 1 for troubleshooting

flagShowRawData = 0;                    % Whether raw (well, raw-ish) data debugging figures need to be shown
flagShowMugshot = 0;                    % Characteristic 
flagShowCorFigure = 1;                  % Raw Correlation and trial-compensated correlation (for every brain)
flagDoTE = 0;                           % Transfer entropy calculation. WARNING: VERY LONG CALCULATION
flagSaveResults = 0;                    % Whether results need to be saved (with set to 0, no risk to overwrite something)

%%% ---------------- Constants ----------------

key = 'cfs';
goodSpikeTimeRange = [0.20 3.5];        % Only this range of times will be kept in spike recordings (to cut out beginning and end artifacts)

nShufflesTE = 1000;                     % How many times TE needs to be reshuffled (1000 was used for the main analysis)

if(flagOneBrainOnly)
    goodBrain = 3;
    iFolder = 1;
    folderName = folderName(3);
end

%%% ---------------- Main cycle ----------------

for(iFolder = 1:length(folderName));
    S = [];     % Empty the data structure
    res = [];   % Empty the output result structure
        
    S = readFiles(folderBaseIn,folderName{iFolder});    % Get all raw data in S.f, and all spike reconstructions in S.s. 1st column sweep#, 2nd column time.
    
    nSweeps = length(S);                                % Number of sweeps
    nCells = size(S(1).dataS,2);                        % Number of cells
    time = S(1).timeS;                                  % Good safe time to use later
    goodSpikeTime = find((time>goodSpikeTimeRange(1)) & (time<=goodSpikeTimeRange(2)));
    lgst = length(goodSpikeTime);
    time = time(goodSpikeTime);                         % Only leave middle of the slides, far enough from artifacts.
    xy = S(1).xy;
    fprintf('%5d cells, %5d stimuli, %5d time points\n',nCells,nSweeps,length(goodSpikeTime));
    if(size(xy,1)>nCells); xy = xy(1:nCells,:); fprintf('Warning: rows(xy)>nCells in data\n'); end;
        
    res.name = folderName{iFolder};         % Copy all important info to the output structure
    res.nCells = nCells;
    
    spikes = [];                            % All sweeps spikes, as a 2D matrix, time down fast, sweeps down slow, cells to the right
    samps = [];                             % Cumulative amplitudes of each response
    fAvC = [];                              % Averages of fluorescence, across cells
    sAvC = [];                              % Averages of spiking, across cells
    data = [];                              % All sweeps spikes, as a 2D matrix, time down, cells right fast, sweeps right slow
    sweepType = mod((1:nSweeps)-1,3)+1;     % Sweep types: 1=C, 2=F, 3=S
    averageShapePerCell = zeros(lgst*3,nCells);  % To keep average traces for every cell, all 3 stim after one another
    for(iSweep=1:nSweeps)
        spikes = [spikes; S(iSweep).dataS(goodSpikeTime,:)];        % Cutting spikes to size and concatenating them in one super-long trace for each cell, for PCAing       
        averageShapePerCell((1:lgst)+(sweepType(iSweep)-1)*lgst,:) = ...
            averageShapePerCell((1:lgst)+(sweepType(iSweep)-1)*lgst,:) + S(iSweep).dataS(goodSpikeTime,:);
        samps = [samps; mean(S(iSweep).dataS(goodSpikeTime,:))];    % Total amplitude of reach cell and each response
        fAvC = [fAvC mean(S(iSweep).dataF(goodSpikeTime,:),2)];     % Average fluorescence across entire OT (all cells lumped) for each stimulus
        sAvC = [sAvC mean(S(iSweep).dataS(goodSpikeTime, :),2)];    % Same for spikes
        data = [data S(iSweep).dataS(goodSpikeTime,:)];             % Unlike in "spikes", concatenation goes right, not down        
    end
    fAvC = bsxfun(@plus,fAvC,-fAvC(1,:));                           % Zero starting points of fluorescence traces
    averageShapePerCell = averageShapePerCell/nSweeps;
    
    %%% Now calculate cross-cell correlations adjusted for averages
    corRaw = corr(spikes(1:end-1,:),spikes(2:end,:))'/nSweeps;      % Flip from normal notation to neuro-notation
    corAdj = reshuffle_corr(spikes,nSweeps)';
    if(flagShowCorFigure)   % Troubleshooting figure for correlations
        figure; 
        subplot(1,3,1); myplot(corRaw); title('raw');
        subplot(1,3,2); myplot(corr(averageShapePerCell(1:end-1,:),averageShapePerCell(2:end,:))); title('on ave');
        subplot(1,3,3); myplot(corAdj); title('adj');
        drawnow();
    end
    
    res.ampF = mean(S(1).dataF(goodSpikeTime,:),1);                 % Typical glowing (fluorescence) of each cell
    traceS = repmat(1:3,1,nSweeps/3)';                  % For each sweep, which stimulus type is it?
    
    % --------------------- Normalize response amplitudes within each triad of stimuli ------------
    sampsN = zeros(size(samps));                        % Normalize amplitude responses within each triad
    for(iTriad=1:floor(nSweeps/3))
        ind = (1:3) + (iTriad-1)*3;        
        sampsN(ind,:) = samps(ind,:)/mean(mean(samps(ind,:)));        
        k = mean(mean(data(:,ind)));                    % Normalization coefficient
        data(:,ind) = data(:,ind)/k;
    end
    
    if(flagShowRawData) % -------------------------- Show raw-ish data of several kinds
        figure('Color','white');    % --- All cells, averaged over time, per response time
        nsp = ceil(sqrt(nCells)/1.5);
        msp = ceil(nCells/nsp);
        for(iCell=1:nCells)
            h = axes('units','norm','Position',[floor((iCell-1)/msp)/nsp mod(iCell-1,msp)/msp 1/nsp 1/msp]); % Instead of subplot(nsp,msp,iCell), which is extremely slow
            tempc = []; temps = []; tempf = [];
            for(iTriad=1:floor(nSweeps/3))
                tempc = [tempc S((iTriad-1)*3+1).dataS(goodSpikeTime,iCell)];
                tempf = [tempf S((iTriad-1)*3+2).dataS(goodSpikeTime,iCell)];
            end
            set(h,'XTickLabel',[],'YTickLabel',[],'nextplot','add');              
            plot( time, mean(tempc,2) , 'b-' , 'parent',h);           
            plot( time, mean(tempf,2) , 'r-' , 'parent',h);
        end                 
        drawnow();
        
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
        
%         figure('Color','white'); hold on; % --- The most raw data view that is barely visible
%         myStyles = 'brg';
%         for(iSweep=1:nSweeps)
%             for(iCell=1:nCells)
%                 plot(time + (iSweep-1)*4, S(iSweep).dataS(goodSpikeTime,iCell)+ (iCell-1)*0.1 , myStyles(mod(iSweep-1,3)+1));
%             end 
%         end
%         title(['All traces; brain ' folderName{iFolder}]);
%         drawnow();
    end    
    
    if(flagShowMugshot)
        figure('Color','white'); % --------------------- mugshot figure starts here ------------
        set(gcf,'Position',get(gcf,'Position')*[1 0 0 0 ; 0 1 0 -0.3 ; 0 0 1 0 ; 0 0 0 1.2]');
        warning('off','all');   % Because scatter ruins everything in Matlab 2014; may get fixed after next update
        
        subplot(3,4,1); % --- Fluorescence
        set(gca,'FontSize',8); hold on;
        plot(time,mean(fAvC(:,sweepType==1),2),'b-'); % C
        plot(time,mean(fAvC(:,sweepType==2),2),'r-'); % F
        plot(time,mean(fAvC(:,sweepType==3),2),'g-'); % S
        title(['Date: ' folderName{iFolder}]);        

        subplot(3,4,2); % --- Spiking
        set(gca,'FontSize',8); hold on;
        plot(time,mean(sAvC(:,sweepType==1),2),'b-'); 
        plot(time,mean(sAvC(:,sweepType==2),2),'r-'); 
        plot(time,mean(sAvC(:,sweepType==3),2),'g-'); 
        title('Spiking');
        legend({'c','f','s'}); legend('boxoff');

        subplot(3,4,9); set(gca,'FontSize',8); % --- Intensity of spiking in each cell        
        scatter(xy(:,1),xy(:,2),30,mean(spikes),'filled');
        xlim([min(xy(:,1)) max(xy(:,1))]);
        ylim([min(xy(:,2)) max(xy(:,2))]);
        title('Intensity');

        subplot(3,4,5); set(gca,'FontSize',8); % --- PCA eigenvalues
        [coeffs,scores,eigs] = princomp(spikes);    % Scores are signals: time down, number right. Coeffs are impacts in each cell (cell down, number right)
        plot(eigs(1:10),'.-');
        title('PCA eigenvalues');

        subplot(3,4,6); set(gca,'FontSize',8); % --- PCA components
        hold on;
        plot(mean(reshape(scores(:,2)+scores(:,1)*0.3,length(time)*3,[]),2),'b-');
        plot(mean(reshape(scores(:,3)+scores(:,1)*0.3,length(time)*3,[]),2),'r-');
        title('Guess at response shapes');

        subplot(3,4,10); set(gca,'FontSize',8); % --- Impacts of two compoments
        rgb = [abs(coeffs(:,3))/max(abs(coeffs(:,3))) zeros(size(coeffs(:,1))) abs(coeffs(:,2))/max(abs(coeffs(:,2)))];
        rgb = 1-repmat(max(rgb,[],2),1,3) + rgb; % Inverses the scale not changing the colors (hard to tell it if it really works)
        scatter(xy(:,1),xy(:,2),30,rgb,'filled');
        xlim([min(xy(:,1)) max(xy(:,1))]);
        ylim([min(xy(:,2)) max(xy(:,2))]);
        title('PCA components 2 and 3');

        subplot(3,4,3); set(gca,'FontSize',8); % --- History of spiking
        myplot(samps); colormap('default');
        title('Amplitude history');
        
        subplot(3,4,4); set(gca,'FontSize',8); % --- History of spiking, adjusted
        [~,ind1] = sort(sweepType);                          % To rearrange the amplitudes by stim type
        [~,ind2] = sort(sum(sampsN,1));
        myplot(sampsN(ind1,ind2)); colormap('default');
        title('Adjusted and sorted');

        subplot(3,4,11); set(gca,'FontSize',8); % --- CS Selectivity
        rgb = [];
        res.ampS = zeros(nCells,3);                                     % Average full responses to CFS will be stored here
        for(iCell=1:nCells)
            thisc = (sampsN(mod(traceS,3)==1,iCell));        
            thisf = (sampsN(mod(traceS,3)==2,iCell));
            thiss = (sampsN(mod(traceS,3)==0,iCell));
            res.ampS(iCell,:) = [mean(thisc) mean(thisf) mean(thiss)];
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
        warning('on','all');
    end
    
    %     doTrent = 0;
    %     if(doTrent) % ----------------------- TRENT tool calculation (doesn't work at all, just freezes the system forever)
    %         caimaging_trent_wrapper(time,data,nCells,traceS);
    %     end
    
    sel = mean(sampsN(mod(traceS,3)==2,:)-sampsN(mod(traceS,3)==1,:),1);  % Crash minus flash, for each cell
    
    res.sel = sel;
    res.xy = xy;
    res.spikiness = mean(spikes);
    res.totalResponses = sampsN;
    res.corRaw = corRaw;
    res.corAdj = corAdj;


    if(flagDoTE)    % ----------------------- Transfer Entropy and other guesses about connectivity       
        
        % -- Simple correlation
        %wc = corr(spikes(1:end-1,:),spikes(2:end,:))';                     % Adjacency needs to be flipped as w is flipped (w_12 is 2 to 1 projection)
        
        % -- Correlation on real data compared to that on surrogate (bootstrapped) data
        sweepMask = reshape(repmat(sweepType(:)',nCells,1),[],1);           % Here sweepType code (1, 2 or 3) is repeated nCell times        
        % wcC = reshuffle_corr(data(:,sweepMask==1),nCells)';                 % Only on crashes
        % wcF = reshuffle_corr(data(:,sweepMask==2),nCells)';                 % Only on scrambles
        % figure; loglog(wcC(:),wcF(:),'.'); xlabel('p-val on C'); ylabel('p-val on F');
        % figure('Color','white'); hist(wcC(:),100); title('Distribution of weights');
        % wct = 1*(wcC>quantile(wcC(:),1-0.01));
        % selectivity_graph(wct,sel,xy(:,1),xy(:,2));
        
        % -- Brute-force binned transfer enthropy            
        teC = transfer_entropy(data(:,sweepMask==1),nCells,3,nShufflesTE); % Loom
        teF = transfer_entropy(data(:,sweepMask==2),nCells,3,nShufflesTE); % Flash
        teS = transfer_entropy(data(:,sweepMask==3),nCells,3,nShufflesTE); % Scrambled
        % figure; plot(wC(:),wF(:),'.'); xlabel('TE on C'); ylabel('TE on F');
        % figure; plot(wC(:),wS(:),'.'); xlabel('TE on C'); ylabel('TE on S');
        
        % figure('Color','white'); hist(w(:),100); title('Distribution of weights');
        % wt = 1*(teC.te>quantile(teC.te(:),1-0.01)); % Look at top 5% of wt
        % selectivity_graph(wt,sel,xy(:,1),xy(:,2));

        res.teC = teC;
        res.teF = teF;
        res.teS = teS;
    end % do TE
    
    fullFileNameOut = [folderBaseOut sprintf('results-%s.mat',folderName{iFolder})];
    fprintf('Accessing output file: %s\n',fullFileNameOut);
    if(exist(fullFileNameOut,'file'))
        fprintf('Uh-oh, this file exists already.\n');
        oldRes = load(fullFileNameOut);
        oldRes = oldRes.res;
        if(flagDoTE)
            fprintf('TE was recalculated, so rewriting the file\n');            
        else
            if(~isfield(oldRes,'teC'))
                fprintf('Old file does not have TE, so rewriting the file\n');                
            else
                fprintf('TE was NOT recalculated, so saving old TE, but updating everything else.\n');
                res.teC = oldRes.teC;
                res.teF = oldRes.teF;
                res.teS = oldRes.teS;
                res.sel = oldRes.sel;
                
%                 res.teC.rawte = res.teC.rawte';   % This section was only ran once on Sep 16 when I moved from slipped TEs (neuro-style) to straight TEs
%                 res.teC.te = res.teC.te';
%                 res.teC.backgroundte = res.teC.backgroundte';
%                 res.teC.p = res.teC.p';
%                 res.teC.corr = res.teC.corr';
%                 
%                 res.teF.rawte = res.teF.rawte';
%                 res.teF.te = res.teF.te';
%                 res.teF.backgroundte = res.teF.backgroundte';
%                 res.teF.p = res.teF.p';
%                 res.teF.corr = res.teF.corr';
%                 
%                 res.teS.rawte = res.teS.rawte';
%                 res.teS.te = res.teS.te';
%                 res.teS.backgroundte = res.teS.backgroundte';
%                 res.teS.p = res.teS.p';
%                 res.teS.corr = res.teS.corr';
            end
        end
    else        
        fprintf('Writing the file.\n');
    end
    if(flagSaveResults)
        save(fullFileNameOut,'res');%,'teF','teS');
    else
        fprintf('Saving flag is set to 0, so not saving the results\n');
    end
    
    if(nargin>0); close all; end;   % If called externally, probably figures aren't needed
end % Each file

end


function Sbag = readFiles(folderBaseIn,folderName)
%%% -------------- Preprocessing -----------------
if(folderName(end)~='\')
    folderName = [folderName '\'];
end
fullFolderName = [folderBaseIn folderName];
fprintf('\nFolder: %s\n',fullFolderName);

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
    fprintf('Reading    mat-file: %s\n',matFileName');
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