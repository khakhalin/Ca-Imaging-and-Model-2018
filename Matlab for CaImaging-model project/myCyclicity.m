function c = myCyclicity(w)
% c= c myCyclicity(w)
% Estimation of cyclicity

if(nargin<1)
    flagTest = 1;
    n = 10;
    w = [ones(n,1) zeros(n,n-1)];       % Let's build a funny random graph with n edges
    w(:) = w(randperm(n^2));
    w = rand(n).*w;
    w = w-diag(diag(w));                % Zeros on the diagonal
    
    w = [0 1 0; 0 0 1; 1 1 0]; n = size(w,1); % For explicit structral testing
else
    flagTest = 0;
end

if(size(w,1)~=size(w,2)); error('Input matrix shoudl be square'); end;

d = 0.9;                                    % dampening factor
maxIterations = 100;                        % d^maxIterations = 3e-5, which feels proper
n = size(w,1);
wp = (ones(size(w))-eye(size(w)))/(n-1);    % Full graph without self-connections of the same size, with weights normalized to sum(in)=1

n = size(w,1);
s = zeros(n);                           % Flow matrix
sp = zeros(n);                          % Flow for the full graph
w = w/max(w(:))/(n-1);                  % Normalize max theoretically possible input current to 1.
%w = w/max(sum(w,1));                   % Another option: normalize max actually possible input current for this graph to 1
for(q=1:maxIterations)
    s  = w' *max(s, eye(n))*d;          % Let the fluid flow; continue pumping 1s in each node, one by one (columns of s)
    sp = wp'*max(sp,eye(n))*d;          % full graph flow
end
c = trace(s)/trace(sp);                 % Ratio stable self-in-flow in our graph compared to a full graph

if(flagTest)
    figure('Color','white');
    G = digraph(w);
    p = plot(G);
    p.NodeCData = diag(s);
    title(c);
end
        
end