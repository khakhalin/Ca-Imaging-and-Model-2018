function c = myCentrality(w,type)
% c = myCentrality(A,type)
%
% Calculates centrality properly
% Expects an adjacency matrix, not a neuro-weight matrix (A = W').

%type = validatestring(type, {'pagerank','reversepagerank'});

n = size(w,1);

switch(type)
    case {'netrank','revnetrank'}               % Modified pagerank, which seems to be identical to Katz centrality
        if(strcmp(type,'revnetrank'))
            w = w';                             % If you want reversed
        end
        d = 0.9;                                % Dampening factor, and a measure of rain (constant flow) to the system (1-d)
        maxIterations = 1000;                   % Max number of iterations
        tol = 1e-4;                             % Tolerance
        cnew = ones(n,1);                       % Starting point
        k = 1/max(sum(w,2));                    % Strongest total output from a node would be set to 1 (some kind of weak scaling)
        for(q=1:maxIterations)
            c = cnew;
            cnew = k*w'*c*d + (1-d);
            if norm(c - cnew, inf) <= tol; break; end;
        end
    case 'pagerank'
        G = digraph(w);
        c = centrality(G,'pagerank');           % Default Matlab pagerank, for troubleshooting
    case 'revpagerank'
        G = digraph(w');
        c = centrality(G,'pagerank');           % Default Matlab pagerank, for troubleshooting
    case 'gatherer'
        c = (1+sum(w,1)')./(1+sum(w,2));        % Sum of inputs divided by sum of outputs (sorta)
    case 'reach'
        reach = breadthdist(w');                % From brain connectivity toolbox: https://sites.google.com/site/bctnet/measures/list
        c = sum(reach,2)/n;                     % Number of nodes that reach to this node (so like Katz but without attenuation)
    case 'revreach'                             % Number of nodes that can be reached from this node (so like reverese Katz)
        reach = breadthdist(w);
        c = sum(reach,2)/n;
    case 'clustering'                           % Clustering coefficient for each node
        c = clustering_coef_wd(w');             
    case {'percol','revpercol'}                 % Some sort of homebrew centrality, as a weighted generalization of reach. I think it's actually also Katz centrality
        if(strcmp(type,'revpercol'))            % https://en.wikipedia.org/wiki/Katz_centrality
            w = w';                             % Wait, I'm silly, it's exactly the same as my netrank, just reversed, and badly calculated. DONT USE IT.
        end                                     % Max possible flow from a point obviously == stable in flow on a reversed graph!
        d = 0.9;                                % dampening factor
        maxIterations = 100;                    % d^maxIterations = 3e-5, which feels proper
        state = zeros(n);
        w = w/max(w(:))/(n-1);                  % Normalize max theoretically possible input current to 1.
        %w = w/max(sum(w,1));                   % Another option: normalize max actually possible input current for this graph to 1
        for(q=1:maxIterations)
            state = w'*max(state,eye(n))*d;     % Let the fluid flow
        end
        c = sum(state,1)';
    case {'assii','assoo','assio','assoi'}      % Assortativities
        ins = sum(w,1)';                        % We assume it is a graph matrix, not the calcculating matrix (graph = calculation'), so w(1,2) is from 1 to 2
        ous = sum(w,2);
        [i,j] = find(ones(size(w)));            % Trick stolen from Sprons' script: now i = 123..n123..n... and j = 111..1222..2...
        switch(type)
            case 'assii'; a = ins(i); b = ins(j);
            case 'assoo'; a = ous(i); b = ous(j);
            case 'assio'; a = ins(i); b = ous(j);
            case 'assoi'; a = ous(i); b = ins(j);
        end
        c = weightedcorr(a,b,w);
end
