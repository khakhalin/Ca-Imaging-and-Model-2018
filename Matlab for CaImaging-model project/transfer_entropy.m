function out = transfer_entropy(data,nCells,nBins,nShuffle)
% te = transfer_entropy(data)
% te = transfer_entropy(data,nBins)
% transfer_entropy(data,nCells,nBins,nShuffle)
%
% Note that all matrices are returned un-flipped (w_ij is i to j), even tho neuro-notation is often backwards (w_ij is j to i).
% In othwer words, G(w) would be a correct graph, but if you want to propagate on this matrix, you'll have to do s = w'*s.
%
% See also: reshuffle_corr, caimaging_structure, selectivity_graph

% My attempt at calculating TE.
% Aug 18 2017: Created
% Aug 22 2017: Reshuffling.
% Aug 23 2017: Better adjustment, p-values, better tester.
% Sep 06 2017: Tried whether linear (not quantiled) binning works better, but it doesn't seem to.
% Sep 16 2017: Reverted from flipped to normal output


if(nargin<1); [data,realw] = fake_data(10); nCells = 10; end;       % Fake data for troubleshooting and the actual w it's based on
if(nargin<3); nBins = 3; end;                                       % Default value for the number of bins
if(nargin<4); nShuffle = 100; end;                                  % Number of surrogate data sets. 10000 is recommended in the lit, but it's super slow

showFigures = 1;            % Whether debugging figures are to be shown
useQuantiles = 1;           % Whether quantile-binning or linear binning is to be used. Recommended: 1

nShuffle = nShuffle+1;                                      % First run for real TEs, so +1
[nTime,rightSize] = size(data);
nSweeps = rightSize/nCells;

% --- Bin each trace into equally-used bins
db = ones(size(data));                      % Binned data
for(ci=1:nCells*nSweeps)
    if(useQuantiles)
        bins = quantile(data(:,ci),(1:nBins-1)/nBins);    
    else
        bins = (1:nBins-1)/nBins*max(data(:,ci));
    end
    for(q=1:nBins-1)        
        db(data(:,ci)>bins(q),ci) = q+1;
    end
end

% --- Now go and brute-force calculate TE for each pair
% We generate nShuffle sets of fake data and then calculate TE for every pair of cells on this fake set.
% Numerically however it's better to calculate some things only once.
te = zeros(nCells,nCells,nShuffle);
for(h=1:nShuffle)
    fprintf('Shuffle %3d/%3d',h,nShuffle);
    d = zeros(nTime*nSweeps,nCells);                        % Reshape data for TE calculation. We'll ignore seams between sweeps to simplify the code
    if(h==1)                                                % If h==1 don't reshuffle the data, but calculate actual TEs        
        for(s=1:nSweeps); d((1:nTime)+(s-1)*nTime,:) = db(:,(1:nCells)+(s-1)*nCells); end;          % Reshape
        dRef = d;                                           % Save this unshuffled set of data for future calculationos
        correlationReference = corr(d(1:end-1,:),d(2:end,:));   % Linear correlation between signals; just to have a reference during visual troubleshooting
    else                                                    % We'll use random permutation of sweeps to reshuffle the data
        badPerm = 1;                                        % But only those permutations work where all sweeps are wrong. So let's assume it's bad...
        while(badPerm)
            perm = randperm(nSweeps);                       % ... and draw new permutations ...                 
            badPerm = sum(perm==(1:nSweeps));               % ... until we get all sweeps paired wrong and badPerm sets to 0
        end        
        for(s=1:nSweeps); d((1:nTime)+(perm(s)-1)*nTime,:) = db(:,(1:nCells)+(s-1)*nCells); end;    % Reshape with reshuffling (sweep s goes to perm(s))
    end
    
    for(j=1:nCells)                         % Cell of interest
        for(i=1:nCells)                     % Potential input to the cell of interest
            if(i==j); continue; end;
            p = zeros(nBins,nBins,nBins);                   % 3D array of binned stats, respectively: this cell time t, this cell t-1, input cells t-1            
            for(t=2:nTime)                                  % For every time tick, populate the 3D histogram
                p(dRef(t,j),dRef(t-1,j),d(t-1,i)) = p(dRef(t,j),dRef(t-1,j),d(t-1,i))+1;  % Here y' (new y) and y (old y) come from actual data, but x comes from shuffled data
            end            
            p = p/sum(p(:));                                % From counts to probabilities: p(y'yx)

            pj_ji = bsxfun(@times,p,1./sum(p,1));           % p(y'|yx)
            pj_ji(isnan(pj_ji)) = 0;                        % replace NANs with zeros, as we'll check for zeros later            
            if(h==1 && (i==1 || (j==1 && i==2)))               % These probabilities need to only be calculated once for each output cell (j)
                pjj = sum(p,3);                                                     % p(y'y)
                pj_j{j}  = repmat(bsxfun(@times,pjj,1./sum(pjj,1)),1,1,nBins);      % p(y'|y) repeated along 3d dimension (x) for simpler calculation later
                pj_j{j}(isnan(pj_j{j})) = 0;        
            end

            for(q=1:nBins^3)
                if(pj_ji(q)==0 || pj_j{j}(q)==0)
                    % When p(y'|yx)==0 we know that p(y'yx) also ==0, and 0*log(0) is assumed to be 0
                    % Same when p(y'|y)==0, 0*log(inf) is assumed to be 0
                else                    
                    te(i,j,h) = te(i,j,h) + p(q)*log(pj_ji(q)/pj_j{j}(q));
                end
            end
        end % i (potential input)
    end % j (cell of interest)
    fprintf('\n');
end % Shuffle

out.rawte = te(:,:,1);                                         
out.te = (te(:,:,1)-mean(te(:,:,2:end),3));                    % Subtract average reshuffled TE from actual TE
out.backgroundte = mean(te(:,:,2:end),3);
for(i=1:nCells)
    for(j=1:nCells)
        out.p(i,j) = max(1,sum(te(i,j,:)>te(i,j,1)))/nShuffle;
    end
end
out.corr = correlationReference;

if(showFigures)
    corSfl = reshuffle_corr(d,nSweeps);         % Shuffled correlation (yet another custom function)
    
    figure('Color','white');
    subplot(2,3,1); myplot(realw); title('Actual connections'); xlabel('Receptient'); ylabel('Sender'); 
    subplot(2,3,2); myplot(correlationReference); title('Correlation'); xlabel('Recepient'); ylabel('Sender');
    subplot(2,3,3); myplot(corSfl); title('Shuffled cor');
    subplot(2,3,4); myplot(te(:,:,1)); title('Raw TE');
    subplot(2,3,5); myplot(out.te); title('Adjusted TE'); xlabel('Recepient'); ylabel('Sender');
    %subplot(2,2,4); myplot(te(:,:,2)); title('Sample wrong TE');
    %subplot(2,3,4); myplot(mean(te(:,:,2:end),3)); title('Average wrong TE');
    %subplot(2,3,5); myplot(-log(out.p)); title('-log(p-value)');
    subplot(2,3,6); myplot(out.te.*(out.p<0.05)); title('TE masked by p<0.05');
    %subplot(2,2,4); myplot(realw);
    
    if(0 & exist('realw','var'))   % If we are in a troubleshooting mode working with generated data
        %figure;
        %subplot(1,2,1); hist(squeeze(te(3,2,:))); hold on; stem(te(3,2,1),1); hold off; title('Existing connection');
        %subplot(1,2,2); hist(squeeze(te(6,7,:))); hold on; stem(te(6,7,1),1); hold off; title('False connection');
        
        figure;
        subplot(2,2,1); plot(realw(:),correlationReference(:),'.'); ylabel('correlation');
        temp = te(:,:,1);
        subplot(2,2,2); plot(realw(:),temp(:),'.'); ylabel('TE');        
        subplot(2,2,3); plot(realw(:),out.te(:),'.'); ylabel('Adjusted TE');        
        %subplot(2,2,4); plot(realw(:),-log(out.p(:)),'.'); ylabel('-log(p-value)');        
        subplot(2,2,4); plot(realw(:),corSfl(:),'.'); ylabel('shuffled cor');        
    end    
end

out.p = out.p;                 % I used to flip all matrices here (with a '), but I do it no more. This is left as a reminder.
out.corr = out.corr;
out.rawte = out.rawte;
out.te = out.te;
out.backgroundte = out.backgroundte;

end


function [data,realw] = fake_data(nCells)

nTime = 1000;
nSweeps = 20;
noiseLevel = 1;

w = zeros(nCells);
% w = abs(randn(nCells)); w = w/max(w(:))*0.9;
% w = rand(nCells)*0.9;
p = primes(numel(w)); for(q=1:length(p)); w(p(q))=1/p(q); end; % Make some of the elements non-zero, but in a predictable way
w(7,3) = 1; % From 3 to 2
w(4,5) = 1; % From 5 to 4
drive = min(10,sum(abs(randn(nCells,3)),2))'/10;     % Shared synaptic drive each cell receives
    
data = zeros(nTime,nCells*nSweeps);
for(iSweep=1:nSweeps)
    cols = (1:nCells)+(iSweep-1)*nCells;
    data(1,cols) = abs(randn(1,nCells))*noiseLevel;     % Randomly initialize
    for(t = 2:nTime)
        temp = (w*data(t-1,cols)')';                    % Propagate, as if spike-vector was a colum-vector (need to rotate as in data time goes down, not cellid)
        data(t,cols) = tanh(1*(temp-median(temp))) + abs(randn(1,nCells))*noiseLevel + sin(t)*drive;    % Inhibition + sigmoid + noise
    end
    data(:,cols) = bsxfun(@plus,data(:,cols),max(0,sin(iSweep*2*pi/nSweeps+(1:nTime)'/nTime*10)-0.5));  % Shared input or artifact
end
realw = w';                                             % Flip from operator w into adjacency matrix w, so that G(realw) would be a correct graph.

% figure; subplot(1,2,1); plot(data(:,1:nCells)); subplot(1,2,2); plot(data(:,(1:nCells)+nCells*floor(nSweeps/2))); % Visualize data

end