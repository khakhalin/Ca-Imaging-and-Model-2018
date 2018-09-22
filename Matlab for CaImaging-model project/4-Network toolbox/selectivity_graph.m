function selectivity_graph(w,sel,trueX,trueY)
% selectivity_graph(w,sel,trueX,trueY)
%
% Vizualizes graph given by the adjacency matrix W using superimposed node properties data SEL.
%
% Note that it expects W to come as an flipped matrix (w_ji goes from i to j), and so promptly flips it,
% to move to standard graph notation where it is w_ij that goes from i to j.
%
% Also still contains the code that compares SEL to various network centrality functions, but it now commented out.

% Aug 20 2017: created.
% Sep 28 2017: revised.

if(nargin<3)    % If true positions aren't given, just fake something, not to break the code below
    nCells = sqrt(numel(w));
    dim = ceil(sqrt(nCells));
    [trueX,trueY] = meshgrid(1:dim,1:dim);    
    trueX = trueX(1:nCells);
    trueY = trueY(1:nCells);
end

% newSel = min(max(sel,-10),10);          % Contain within [-10 10]
% newSel = newSel(:);
newSel = sel;

figure('Color','white'); % --- Figure showing graph topology
w = floor(w'+0.5);              % Invert as w_ij is a projection from j to i, and then threshold the heck out of it
G = digraph(w);                 % Directed graph

subplot(1,2,1);
pbaspect([1 1 1])
p = plot(G,'NodeLabel',{},'XData',trueX(:),'YData',trueY(:));
p.NodeCData = newSel;
axis('off'); title('Original coordinates');

subplot(1,2,2);
good = (sum(w,1)>0 | sum(w',1)>0);  % Let's not plot unconnected nodes
p = plot(digraph(w(good,good)),'NodeLabel',{});
p.NodeCData = newSel(good);
axis('off'); title('Optimized coordinates');

% % --- Different measures of centrality, and whether they correlate with selectivity
% 
% % From Matlab 'centrality' help (their own script that works for digraphs). Basically except for betweenness it seems to be rather useless.
% % And no code is available, as they are all written as hidden MEXX files or something like that.
% %         'outdegree' - number of successors of node i.
% %          'indegree' - number of predecessors of node i.
% %      'outcloseness' - inverse sum of distances from node i to all reachable nodes. So ==0 for an unconnected node. Ignores weights.
% %       'incloseness' - inverse sum of distances from all nodes to node i, if node i is reachable from these nodes. Same thing; ignores weights.
% %       'betweenness' - number of shortest paths between other nodes that pass through node i. Ignores weights.
% %          'pagerank' - ratio of time spent at node i while randomly traversing the graph. Also ignores weights.
% %              'hubs' - nodes with successors of high authority scores. And both of these two also ignore weights.
% %       'authorities' - nodes with predecessors of high hub scores.
% 
% hF2 = figure('Color','white'); 
% % hF3 = figure('Color','white'); 
% 
% % centralityType = {'indegree', 'outdegree', 'incloseness', 'outcloseness', 'betweenness', 'pagerank','hubs', 'authorities'}; % That's for evil Matlab
% centralityType = {'pagerank','revpagerank','netrank','revnetrank','gatherer','reach','revreach','clustering'};
% for(q=1:length(centralityType))
%     %cent = centrality(G,centralityType{q});    % Evil Matlab centrality
%     cent = myCentrality(w,centralityType{q});  % My nice proper home-grown centrality
%     figure(hF2); subplot(2,4,q);
%     plot(newSel,cent,'.');
%     titleString = centralityType{q};
%     [rho,pval] = corr(newSel,cent);
%     if(pval<0.05)
%         pol = polyfit(newSel,cent,1);
%         hold on;
%         xl = get(gca,'XLim');
%         plot(xl,polyval(pol,xl),'r-');
%         hold off;
%         titleString = sprintf('%s\np=%s, r=%s',titleString,myst(pval),myst(rho));
%     end
%     title(titleString);
%     
%     if(0) % To visually check that centrality measures are calculated properly
%         figure(hF3);
%         subplot(2,4,q);
%         pbaspect([1 1 1])
%         p = plot(G,'NodeLabel',{}); %,'XData',trueX(:),'YData',trueY(:));
%         p.NodeCData = cent;
%         axis('off'); title(titleString);
%     end
% end

end