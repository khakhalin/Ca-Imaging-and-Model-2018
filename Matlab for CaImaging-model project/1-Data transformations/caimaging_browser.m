function caimaging_browser(S,mode,key)
% CAIMAGING_BROWSER(S)
% CAIMAGING_BROWSER(S,mode)
% CAIMAGING_BROWSER(S,mode,key)
%
% A function for browsing through Ca imaging data.
%
% Here MODE can be 'f' for f-data, or 's' for spike data.
% KEY provides a key for the stimuli (a row of several characters, like 'cfs' for example; 
% stimuli are assumed to be circular).

% Aug 22 2014: Key added.

showSpikes = 0;
if(nargin>1)
    if(strcmp(mode,'s'))
        showSpikes = 1;
    end
end
if(nargin<3)
    key = [];
end

nSweeps = length(S);
cellData = [];
for(iSweep=1:nSweeps)
    if(isempty(cellData))                   % First time
        if(showSpikes)
            [n,m] = size(S(iSweep).dataS);
        else
            [n,m] = size(S(iSweep).dataF);
        end
        nCells = m;
        for(iCell=1:nCells)
            cellData{iCell} = [];           % Creating bags
        end
    end    
    for(iCell=1:nCells)        
        if(showSpikes)
            cellData{iCell} = concatenan(cellData{iCell},S(iSweep).dataS(:,iCell));
        else
            cellData{iCell} = concatenan(cellData{iCell},S(iSweep).dataF(:,iCell)/S(iSweep).dataF(1,iCell));
        end
    end
end

if(~isempty(key))                           % If key was provided, need to average the data
    for(iCell=1:nCells)
        [nTime,nTrials] = size(cellData{iCell});
        for(iKey = 1:length(key))
            temp(:,iKey) = mean(cellData{iCell}(:,mod(1:nTrials,length(key))==iKey-1),2);
        end
        cellData{iCell} = temp;
    end
end

scatterbrowser(S(1).xy(:,1), -S(1).xy(:,2), cellData);

end