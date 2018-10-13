function c = myNetMeasure(w,type)
% c = myNetMeasure(A,type)
% 
% Calculates different measures for the entire network (not a value per node, but one value for the entire network)
% Takes adjacency matrix A (a_ij) rather than W (w_ji).

% Oct 13 2018: former myCentrality split into centralities proper and this one.

[n,m] = size(w);

switch(type)
    case 'reciprocity'                          % Returns something like a weighted share of reciprocal connections
        if(min(w(:)<0))                         % Currently doesn't really support negative weights
            warning('Reciprocity: adjacency matrix contains negative weights. It probably makes the outputs meaningless');
        end
        if(n~=m)
            error('Reciprocity needs the matrix to be symmetric');
        end
        w = w.*(ones(n)-eye(n));                % Remove diagonal
        w = w/max(w(:));                        % Normalize
        c = sum((w.*w'))/sum(w);
        
    case {'assii','assoo','assio','assoi'}      % Assortativities (one number, not a vector!)
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
    otherwise
        error('Unknown measure. Check your spelling?');
end

end