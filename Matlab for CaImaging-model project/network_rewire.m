function w = network_rewire(w,m)
% w = network_rewire(w,m)
% w = network_rewire(w)
% Maslov degree-preserving cross-wiring.
% m is the number of shuffles to perform (default m=3*Nedges)

% Oct 02 2017: Verified as functional.
% Aug 06 2018: Updated so that zero-edges are now "swapped", as it doesn't do anything, which made it disfunctional on sparce graphs

if(nargin<1) % Test
    w = eye(100);
    w = [w(end,:); w(1:end-1,:)];
    figure; myplot(network_rewire(w));
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
        if(w(i1,j1)~=0 || w(i2,j2)~=0)                  % No need to "swap" non-connections, as nothing changes when you do that
            direct = w(i2,j2);    w(i2,j2) = w(i2,j1);    w(i2,j1) = direct;    % Cross become directs, while direct become crosses
            direct = w(i1,j1);    w(i1,j1) = w(i1,j2);    w(i1,j2) = direct;
            count = count+1;
        end
    end
end

end