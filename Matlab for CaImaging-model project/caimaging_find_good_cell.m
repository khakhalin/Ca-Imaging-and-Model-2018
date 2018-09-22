function res = caimaging_find_good_cell(S,singleFigure)
% cellid = CAIMAGING_FIND_GOOD_CELL(S)
% cellid = CAIMAGING_FIND_GOOD_CELL(S,suppressNewFigure)
%
% Helps to find a good amplitude cell. Returns a suitable number, and also displays a neat plot.
% suppressNewFigure may be set to 1.

% Apr 29 2014: Created
% Mar 20 2015: Now also returns default "good cell" from the middle of the list.

if(nargin<2)
    singleFigure = 0;
end

depthIntoS = length(S);     % How many trials to average
nShow = 10;                 % How many best trials to show

[nTime,nCells] = size(S(1).dataF);
nTime = nTime-20;                   % To be on a safe side, let's ignore possible variations in length across trials.

avResp = zeros(nTime,nCells);

for(is = 1:depthIntoS)
    avResp = avResp + S(is).dataF(1:nTime,:);       % Average responses    
end

avResp = bsxfun(@plus,avResp,-avResp(1,:));         % Remove biases

[~,ind] = sort(var(avResp,[],1));                   % Find biggest ones

avResp = avResp - min(avResp(:));
avResp = avResp / max(avResp(:));

subset = ind(end-nShow+1:end);                      % Take only topp candidates
res = subset(round(nShow/2));                       % Return one from the middle, it's probably fine.

if(singleFigure~=1)                                 % Figure creation may be suppressed
    figure;
end
plot(avResp(:,subset)+ones(nTime,nShow)*diag(1:nShow)*0.2);
set(gca,'YTick',(1:nShow)*0.2+avResp(1,1),'YTickLabel',{subset});

title('Strongest cells in this set');

end