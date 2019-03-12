function [clustInd,modHistory] = spectralClustering(w,maxClustNumber,sigma,minClustNumber)
% [ind,modHistory] = spectralClustering(weightMatrix,maxClustNumber,sigma,minClustNumber)
% [ind,modHistory] = spectralClustering(weightMatrix,maxClustNumber,sigma)
% [ind,modHistory] = spectralClustering(weightMatrix,maxClustNumber)
% [ind,modHistory] = spectralClustering(weightMatrix)
% ind = spectralClustering(signals)
%
% Based on Ng-Jordan-Weiss spectral clustering algorithm. Tries different numbers of clusters 
% between 1 and MAXCLUSTNUMBER, calculates Newman's spectral modularity for each split, and picks the one with
% the highest modularity. If MAXCLUSTNUMBER is omitted, defaults to 7. If MINCLUSTNUMBER isn't given, defaults to 1.
% SIGMA (default 100) is a parameter for NG-Jordan-Weiss approach, for exponentially transforming 
% Eucledian distances (higher numbers = more attention to weak edges). The main input W is a positive-only correlation matrix.
% 
% Returns a column of cluster indices CLUSTIND, and a history
% of modularity measurements for different cluster numbers (MODHISTORY).
%
% If run without any inputs, tests itself.

% Aug 21 2018: Adopted and adapted

% Reference: Ulrike von Luxburg, "A Tutorial on Spectral Clustering", Statistics and Computing 17 (4), 2007
%
% This is an implementation of:
% Ng, A., Jordan, M., and Weiss, Y. (2002). On spectral clustering: analysis and an algorithm. In T. Dietterich,
% S. Becker, and Z. Ghahramani (Eds.), Advances in Neural Information Processing Systems 14 
% (pp. 849 – 856). MIT Press.
% Adapted from the code of Asad Ali, GIK Institute of Engineering Sciences & Technology, Pakistan; Email: asad_82@yahoo.com
%
% The modularity is based on Newman 2006, and according to this paper, it works for positive-only weighted graphs as well
% (but would require an adjustement for a non-positive matrix): Gómez, S., Jensen, P., & Arenas, A. (2009). Analysis of 
% community structure in networks of correlated data. Physical Review E, 80(1), 016114.


if(nargin<1) % Test
    nComponents = 27;
    nNodes = 41;
    blockSize = 7; % Shoudl be much smaller than both nNodes and nComponents
    %signalTypes = sin(bsxfun(@times,(1:100)'/10,rand(1,nComponents)*10));
    %signalTypes = sin(bsxfun(@times,(1:300)'/10,(1:nComponents)));
    signalTypes = randn(300,nComponents);
    trueWeights = rand(nComponents,nNodes)*0.0+0.9;
    trueWeights = trueWeights.*(ones(size(trueWeights))*0.1 + 0.9*blkdiag(ones(blockSize),ones(blockSize),ones(nComponents-blockSize*2,nNodes-blockSize*2))); % Introduce blocks
    %trueWeights = ones(size(trueWeights))*0.1 + 0.9*blkdiag(ones(blockSize),ones(blockSize),ones(nComponents-blockSize*2,nNodes-blockSize*2)); % Introduce blocks
    %trueWeights = blkdiag(ones(blockSize),ones(blockSize),ones(nComponents-blockSize*2,nNodes-blockSize*2)); % Introduce blocks
    trueWeights = trueWeights(:,randperm(nNodes));          % Reshuffle
    signals = signalTypes*trueWeights;
    signals = signals+randn(size(signals))*0.0;
    corMatrix = corr(signals);
    %corMatrix = 1*(corMatrix>0.5);                          % Threshold. Makes everything easier.
    
    hF = figure; subplot(2,2,1); myplot(trueWeights); subplot(2,2,2); myplot(corMatrix);
    
    [clustInd,modHistory] = spectralClustering(corMatrix,nNodes-1,100);    
    [~,ind] = sort(clustInd);
    
    figure(hF);
    subplot(2,2,3); plot(modHistory,'.-'); title('Modularity(Nclust)');
    subplot(2,2,4); myplot(corMatrix(ind,ind));
    return
end

%%% ------------------ Now the actual function ----------------

if(nargin<4)
    minClustNumber = 1;
end
if(nargin<3)
    sigma = 100;                     % 1 in the original sample; 100 or 1000 seems to work better here
end
if(nargin<2)
    maxClustNumber = 7;
end

if(min(w(:))<0)
    fprintf('Warning: some edges are negative. Switching to absolute values.\n');
    w = abs(w);                         % Make the input matrix all-positive (it should be, but just in case let's make sure it is)
end
nNodes = size(w,1);
maxClustNumber = min(maxClustNumber,nNodes);    % Cannot have more clusters than nodes

% [d(x,y)]^2 = 2-2*cor(x,y)  https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering
affinity = 2*(ones(size(w))-w);                 % From correlations to distances
affinity = exp(-affinity/sigma);                % Required exp transformation

% figure; subplot(1,2,1); imshow(affinity,[]); subplot(1,2,2); plot(affinity);

degreeMat = diag(sum(affinity));

% compute the normalized laplacian / affinity matrix (method 1)
%NL1 = D^(-1/2) .* L .* D^(-1/2);
for i=1:nNodes
    for j=1:nNodes
        NL1(i,j) = affinity(i,j) / (sqrt(degreeMat(i,i)) * sqrt(degreeMat(j,j)));  
    end
end

% perform the eigen value decomposition
[eigVectors,~] = eig(NL1);
if(~isreal(sum(eigVectors)))
    fprintf('WARNING: complex eigenvalues! Everything may be incorrect! (Looking at real parts only)\n');
    eigVectors = real(eigVectors); 
end

m = sum(w(:))/2;                                                % Total number of edges = sum all degrees / 2
degProd = diag(sum(w,1))*ones(nNodes)*diag(sum(w,2))/(2*m);     % a matrix where each element = k^out_i * k^in_j / sum(k), where k_i is out-degree of node i. 
                                                                % and k_j is in-degree of node j. For Newman's modularity
modHistory = zeros(maxClustNumber,1);
for(k=minClustNumber:maxClustNumber)                                            
   [~,modHistory(k),~] = doClustering(eigVectors,w,m,degProd,k);
end
bestK = find(modHistory==max(modHistory));                      % k with highest modularity

[clustInd,~,dist] = doClustering(eigVectors,w,m,degProd,bestK);   % Last (optimal) run

end


function [clustInd,modul,dist] = doClustering(eigVectors,w,m,degProd,nClust)
% Working horse of a subroutine
nNodes = size(w,1);

% select k largest eigen vectors
nEigVec = eigVectors(:,(size(eigVectors,1)-max(1,nClust-1)): size(eigVectors,1));

% construct the normalized matrix U from the obtained eigen vectors
U = [];
for i=1:size(nEigVec,1)
    n = sqrt(sum(nEigVec(i,:).^2));    
    U(i,:) = nEigVec(i,:) ./ n; 
end

[clustInd,~,sumd] = kmeans(U,nClust);
dist = sum(sumd);

% Now  modularity calculation after Newman 2006 (Modularity and community structure in networks)
s = zeros(nNodes);                          % Future index matrix. Zero for those nodes that don't belong to the same block.
for(i=1:nClust)  % For every cluster        
    s(clustInd==i,clustInd==i) = 1;         % 1 for those pairs of nodes that belong to the same block.        
end
modul = sum(s(:).*(w(:)-degProd(:)))/(4*m);         % Newman's formula

end