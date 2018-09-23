function w = network_rewire(w,m)
% w = network_rewire(w,m)
% w = network_rewire(w)
% Degree-preserving cross-wiring.
% (Out-degree only for weighted graph; both out- and in- for binary directed graph).
% W is the weight matrix, understood in operator (flipped) way: W_ji is an edge from i to j.
%
% M is the number of shuffles to perform (default m=3*Nedges)
%
% Based on: Maslov, S., & Sneppen, K. (2002). Specificity and stability in topology of protein networks. Science, 296(5569), 910-913.
% But generalized to weighted networks (so it cannot rely on many weights being ==0).

% Oct 02 2017: Verified as functional.
% Aug 06 2018: Updated so that zero-edges are now "swapped", as it doesn't do anything, which made it disfunctional on sparce graphs
% Sep 23 2018: Different logic to approximate binary graphs.

if(nargin<1) % Test
    fprintf('Testing network rewire. In both pairs, both lines should be the same.\n');
    w = floor(0.1+rand(100));
    fprintf('In-degrees:\n');
    fprintf('%2d ',histcounts(sum(               w,1),(1:40)-0.5)); fprintf('\n');
    fprintf('%2d ',histcounts(sum(network_rewire(w),1),(1:40)-0.5)); fprintf('\n');
    fprintf('Out-degrees:\n');
    fprintf('%2d ',histcounts(sum(               w,2),(1:40)-0.5)); fprintf('\n');
    fprintf('%2d ',histcounts(sum(network_rewire(w),2),(1:40)-0.5)); fprintf('\n');
    w = 0; % (To suppress output)
end

n = size(w,1);
nEdges = sum(w(:)~=0);

if(nargin<2)
   m = 3*nEdges;
end

count = 0;
while(count<m)
    i1 = floor(rand(1)*n)+1;
    j1 = floor(rand(1)*n)+1;
    i2 = floor(rand(1)*n)+1;
    j2 = floor(rand(1)*n)+1;
    if((i1~=j1) && (i2~=j2) && (i1~=j2) && (i2~=j1))    % Don't mess with self-connections in either direction; these should always stay 0
        if(w(j1,i1)>w(j1,i2) && w(j2,i2)>w(j2,i1))                  % No need to "swap" non-connections, as nothing changes when you do that
            direct = w(j2,i2);    w(j2,i2) = w(j2,i1);    w(j2,i1) = direct;    % Cross become directs, while direct become crosses
            direct = w(j1,i1);    w(j1,i1) = w(j1,i2);    w(j1,i2) = direct;
            count = count+1;
        end
    end
end

end