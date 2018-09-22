function A = create_modular_network(nCells, nModules, pIn, pOut)
% A = create_modular_network(nCells, nModules, pIn, pOut)
% A = create_modular_network(nCells, nModules)
% A = create_modular_network(nCells)
% 
% Creates a modular network of NCELLS nodes, with NMODULES (default value ~ sqrt(nCells)).
% The cells are more connected within the modules (pIn, default value of 0.5) than between the modules 
% (pOut, default of 0.2), and while connections within the modules are random, connections between
% the modules are turning them into layers of a hierarcical network. Checks that there are no isolated nodes,
% but doesn't guarantee that the graph is all in one piece (weakly connected).
% 
% If run without parameters, tests itself.

% Aug 16 2018: Created

if(nargin<1)    % Test
    A = create_modular_network(50);
    figure('Color','white'); plot(digraph(A));
    return;
end
if(nargin<2); nModules = floor(sqrt(nCells)); end
if(nargin<3); pIn = 0.6; end
if(nargin<4); pOut = 0.1; end


w = floor(rand(nCells)+pIn);                            % Prototype matrix
w = max(0,w-w');                                        % Remove reciprocal (symmetrical) connections and self-connections
A = zeros(nCells);                                      % Place where the final matrix will be
groups = mod((1:nCells)-1,nModules)+1;                  % Module id for each cell
for(iMod=1:nModules)
    g = (groups==iMod);                                 % Index vector for current group    
    A(g,g) = w(g,g);                                    % Copy these (within), but leave them disconnected
    gNext = (groups==(1+mod(iMod,nModules)));           % Next group, in a cycle
    A(g,gNext) = floor(rand(sum(g),sum(gNext))+pOut);   % Connect this group to the next group
end

isolates = find(sum(A,1)==0 & sum(A,2)'==0);            % Find isolated points (if any)
if(length(isolates)>0)                                  % If ther eare some...
    A(isolates,1+mod(isolates,nCells)) = 1;             %  ... connect them to something. At least they won't be lying there as solitary nodes.
end

end