function [flow,revFlow,flowPerNode] = network_flow(w)
% [flow,revFlow,flowPerNode] = network_flow(w)
%
% Essentially it calculates Katz centrality.
%
% flow: one number, max(KatzRank)-mean(KatzRank)
% revFlow: exactly the same thing, but for w'


% Aug 31 2017: created
% Jul 11 2018: remodeled to a more sane model
% Jul 17 2018: I accepted that it's just Katz centrality. Edge-wise Flow is not reported back for now, as I don't have use for it.

if(nargin<1)
    fprintf('No inputs, so running in test mode.\n'); w = test_matrix(); flagTest = 1;  % If no inputs, run a test
else
    flagTest = 0;
end

n = size(w,1);      % Number of nodes

%%% --- Sharpness
% %sh = max(w,[],2)./sum(w,2);            % Sharpness of the output (max output divided over total output)
% sh = skewness(w,0,2);                   % Skewness, adjusted for bias
% sh(isnan(sh)) = 0;                      % If it got divided by 0 because the node is a dead-end, set it to 0. Or if we got skewness of uniform.
% % figure('Color','white'); plot(digraph(w.*(w>0.7)),'NodeLabel',{},'NodeCData',sh); % A figure to check whether the sharpness measure makes sense
% shi = repmat(sh(:)',n,1);               % Matrix of properties for sending neurons
% sho = shi';                             % Matrix of properties for receiving neurons
% 
% sa = weightedcorr(shi(:),sho(:),w(:));  % Sharpness-based mixing: weighted correlation of sharpness on both ends of each edge

%%% --- Stable flow
[vFlow,eFlow] = calculate_flow(w);        % vFlow - flow into each VERTEX (aka pagerank'); eFlow - flow across each EDGE
flow =    max(vFlow)- mean(vFlow);                      % Global reach via flow (Katz centrality)
flowPerNode = vFlow;                                    % pageranks themselves (if anybody asked)

if(nargout>1 || flagTest)
    [vFlowP,eFlowP] = calculate_flow(w');
    revFlow = max(vFlowP)-mean(vFlowP);                     % Reverse global reach
end

if(flagTest)
    G = digraph(w);
    GP = digraph(w');
    Gfl = digraph(eFlow);
    GflP = digraph(eFlowP);
    %G = digraph(flowMat);
    figure('Color','white'); 
    
    subplot(2,2,1);    
    edgeScale = 3/max(Gfl.Edges.Weight(:));
    plot(G,'NodeLabel',{},'NodeCData',vFlow,'LineWidth',Gfl.Edges.Weight*edgeScale); set(gca,'visible','off');    
    subplot(2,2,2);    hist(vFlow,0:0.01:2);    xlim([0 2]);
    hold on; plot(mean(vFlow),1,'ro');
    
    subplot(2,2,3);    
    edgeScale = 3/max(Gfl.Edges.Weight(:));
    plot(GP,'NodeLabel',{},'NodeCData',vFlowP,'LineWidth',GflP.Edges.Weight*edgeScale); set(gca,'visible','off');    
    subplot(2,2,4);    hist(vFlowP,0:0.01:2);   xlim([0 2]);
    hold on; plot(mean(vFlowP),1,'ro');
    
    drawnow();    
end

end


function [vFlow,eFlow] = calculate_flow(w);

d = 0.9;                                % Dampening factor, and a measure of rain (constant flow) to the system (1-d)
maxIterations = 100;                    % Max number of iterations
n = size(w,1);
vFlow = ones(n,1);                      % Starting point
k = 1/max(sum(w,2));                    % Strongest total output from a node would be set to 1 (some kind of weak scaling)
for(q=1:maxIterations)
    vFlow = k*w'*vFlow*d + (1-d);    
end
eFlow = k*w.*repmat(vFlow,1,n);         % Flow across each edge

end


function w = test_matrix()
% Generates a test matrix.

% This one looks like a Y and a singlish node on the side, in a weak cycle
% w = [0  0  1  0  0  0
%      0  0  1  0  0  0
%      0  0  0  1  0  0
%      0  0  0  0  1  .1
%      0  0  0  0  0  0
%      0  .1  0  0  0  0]; 
 
w = [0  0  1  0  0  0
     0  0  1  0  0  0
     0  0  0  1  0  0
     0  0  0  0  1  0
     0  1  0  0  0  0
     0  0  0  1  0  0]; 
 
 
 % w = w+abs(randn(size(w))*0.2);     % Noisify (if needed)
 w = w.*(1-eye(size(w,1)));         % Remove diagonal

end