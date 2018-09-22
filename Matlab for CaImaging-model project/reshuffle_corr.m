function [w,pval,data2] = reshuffle_corr(data,nSweeps,delay)
% [w,pval,unbiasedData] = reshuffle_corr(data,nSweeps,delay)
% [w,pval,unbiasedData] = reshuffle_corr(data,nSweeps)
%
% Calculates correlations between straight and reshuffled data, and returns
% effect sizes.
%
% DATA: All sweeps, as a 2D matrix, time down fast, sweep i down slow, cells to the right
% NSWEEPS: the number of sweeps. We could have used a 3D matrix, but a series of 2D matrix is easier to imagine
% DELAY: by how many frames to offset the second track during correlation. Defaults to 1,
%   or 1 frame delay, which means that it is looking for CAUSAL correlation, not real correlation.
%   that's for compatibility. Set to 0 for "normal" correlation.
%   Note that for now there's a minor simplification for DELAYs!=0 that adds some trivial amount of noise to the correlation.
%       it shouldn't matter for normal signals, but if your signals are super-short and there are lots of them, yo may want to rewrite this part.
%
% Returns: W adjustsed correlation matrix, PVAL for these correlations, and UNBIASEDDATA for any other processing.

% Aug 21 2017: created
% Sep 16 2017: now returns both p and correlation values.
% Aug 17 2018: updated to "normal', not causal correlation.
% Aug 20 2018: transitioned from shuffling to calculation on adjusted dataset, as with linear correlations shuffling is useless

if(nargin<1) % Testing mode
    myNSweeps = 13;
    backgroundProfile = sin((1:20)/5)';
    hum = repmat(backgroundProfile,myNSweeps,7)*diag(rand(7,1));            % myNSweeps sweeps, 7 cells
    signalProfile = rand(20*myNSweeps,1);                                   % True signal
    trueSignal = repmat(signalProfile,1,7)*diag([1 0 .8 0 0 0 0]);   % So true correlation is only between 1 and 3
    %figure; myplot(trueSignal);
    fullSignal = hum+trueSignal;
    fullSignal = fullSignal + randn(size(fullSignal))*0.1;
    %figure; myplot(fullSignal);
    figure; subplot(2,2,1); plot(hum); subplot(2,2,2); plot(fullSignal);
    subplot(2,2,3); myplot(corr(fullSignal)); title('Normal correlation');
    [c,p] = reshuffle_corr(fullSignal,myNSweeps,0);
    subplot(2,2,4); myplot(c); title('Reshuffled correlation');
    return
end

if(nargin<3); delay = 1; end

[n,nCells] = size(data);
nTime = n/nSweeps;
pval = ones(nCells);                                % pre-alocate; with ones, as they are p-values
w = zeros(nCells);                                  % Actual correlations

data2 = zeros(size(data));                          % Adjusted data
avResponse = zeros(nTime,nCells);                   % Average response

for(iSweep=1:nSweeps)                                           % For each sweep:
    t = (1:nTime) + (iSweep-1)*nTime;                           % find time indices that correspond to iit
    data2(t,:) = bsxfun(@plus,data(t,:),-mean(data(t,:)));      % unbias
    data2(t,:) = bsxfun(@times,data2(t,:),1./std(data2(t,:)));  % normalize
    avResponse = avResponse + data2(t,:);                       % add to the running average
end
avResponse = avResponse/nSweeps;                                % computer average
data2 = data2-repmat(avResponse,nSweeps,1);                     % subtract average from each instance of the the data

if(delay==0)
    [w,pval] = corr(data2);
else
    [w,pval] = corr(data(1:(end-delay),:),data((1+delay):end,:));   
    % A bit of an approximation, as transitions from one sweep to another shoudl technically be taken out.
    % For now I have no motivation to code this more carefully, as I don't plan to use it with delay!=0 in the nearest future, and even if I do,
    % the amount of noise it would add is likely to be minimal.
end

end



function old_function()
% Old abandoned obsolete version. This sort of logic (random permutations) works well for non-linear computations, but it is
% completely unnecessary in this case, as corralations are of course linear.

nShuffle = 50;              % Number of extra reshuffle runs. Recommended: 10
flag_pics = 0;              % Debugging pictures

[n,nCells] = size(data);
nTime = n/nSweeps;
pval = ones(nCells);                                % pre-alocate; with ones, as they are p-values
w = zeros(nCells);                                  % Actual correlations

corRaw  = zeros(nCells);
corFake = zeros(nCells,nCells,nShuffle);
time=1:(nTime-delay);                               % One tick less than full time, as we'll be shifting by DELAY
for(iShuffle=1:nShuffle)                            % How many times to shuffle
    perm = safeperm(nSweeps);                       % Random permutation that shares no digits with a proper 1:n sequence
    for(iSweep=1:nSweeps)                           % We used built-in vector stuff to run all cells at once, so we have to manually cycle through sweeps
        shift1 = (iSweep-1)*nTime;                  % Location of this sweep in the DATA variable
        shift2 = (perm(iSweep)-1)*nTime;            % Wrong of some other random sweep within the DATA variable
        if(iShuffle==1)                             % Proper correlation only needs to be calculated once, on the first run, so it's skipped on all other runs          
            corRaw =  corRaw + corr(data(shift1+time,:),data(shift1+time+delay,:)); % Pair proper x with proper y, calculate corr for all cells in none gow
        end        
        corFake(:,:,iShuffle) = corFake(:,:,iShuffle)+ corr(data(shift1+time,:),data(shift2+time+delay,:));    % Pair proper x with wrong y
    end    
end
corRaw = corRaw/nSweeps;                            % From a sum to an average
corFake = corFake/nSweeps;

for(i=1:nCells)
    for(j=1:nCells)
        [~,pval(i,j)] = ttest(squeeze(corFake(i,j,:))-corRaw(i,j));
    end
end

w = corRaw-mean(corFake,3);                         % mean across the 3d dimension

end


function p = safeperm(n)
% Returns a random permutation that shares no digits with proper sequence. Lazy programming (just sample randperm until we get a good one).
badPerm = 1;                                        % But only those permutations work where all sweeps are wrong. So let's assume it's bad...
while(badPerm)
    p = randperm(n);                       % ... and draw new permutations ...                 
    badPerm = sum(p==(1:n));               % ... until we get all sweeps paired wrong and badPerm sets to 0
end
end



function obsolete_function()
data = bsxfun(@plus,data,-mean(data,1));            % Unbias
data = bsxfun(@times,data,1./std(data,[],1));       % Normalize
data = bsxfun(@plus,data,-mean(data,1));            % ReUnbias

% -- Testing section; just making sure the data came in right
% figure; hold on;
% cl = 'brkgmc';
% for(iSweep=1:nSweeps); plot(data(:,(iSweep-1)*nCells+(1:nCells))+iSweep,cl(mod(iSweep-1,length(cl))+1)); end;
% for(iCell=1:nCells); plot(data(:,(0:(nSweeps-1))*nCells+iCell)+iCell,cl(mod(iCell-1,length(cl))+1)); end;
% error();

s = zeros(nShuffle,1);                  % To keep reshuffled correlations
for(i=1:nCells)    
    fprintf('%3d/%3d ',i,nCells);    
    for(j=1:nCells)
        if(i==j); continue; end                 % No self-connections so just keep w_ii==0
        c = corr( reshape(data(1:(end-1),i+(0:(nSweeps-1))*nCells),[],1) , ...
                  reshape(data(2:end,    j+(0:(nSweeps-1))*nCells),[],1) );     % Calculate i to j transfer correlation        
        for(iShuffle=1:nShuffle)
            a = data(1:(end-1),i +         (0:(nSweeps-1))*nCells);
            b = data(2:end,    j + shuffle((0:(nSweeps-1))*nCells));            
            %s(iShuffle) = corr(a(:),b(:));
            s(iShuffle) = a(:)'*b(:)/(numel(a)-1);            
        end
        p = max(sum(s>c)/length(s),1/nShuffle);                       % Probability estimation
        if(flag_pics)
            figure; hold on; hist(s,100); stem(c,1,'r'); hold off; set(gca,'FontSize',6); title(myst(p)); drawnow; % error();            
        end
        %fprintf(' %s',myst(sum(s>c)/length(s)));
        pval(i,j) = p;
        w(i,j) = c-mean(s);        
        if(0)   % Report significance on the fly
            if(p<0.001); fprintf('i'); elseif(p<0.01); fprintf(':'); elseif(p<0.05); fprintf('.'); else; fprintf('_'); end;
            if(mod(j,80)==0); fprintf('\n        '); end;
        end
    end
    fprintf('\n');    
end

%w(:) = -log(w(:));
%w(:) = -log(w(:)).*fdr(w(:));                   % log-inverse, but zero those that aren't significant

end


function s = shuffle(a) % Only works for vectors
s = a(randperm(length(a)));
end
