function s = myst(v,shortflag)
% S = MYST(V)
% S = MYST(V,shortflag)
%
% MYST is a short for My String. It creates a nice intuitive p-value string for V p-value.
% (fixed-point for large numbers, exponential for small ones)
% If SHORTFLAG==1, it uses . .. * ** coding instead.

% Nov 05 2012: review
% Apr 25 2015: short form.
% Jul 02 2015: minor adjustment.

if(nargin<2)
    shortflag = 0;
end

if(~shortflag)
    if(abs(v)<0.0001)
        s = sprintf('%1.0e',v);
    elseif(abs(v)<0.01)
        s = sprintf('%1.4f',v);
    elseif(abs(v)<0.1)
        s = sprintf('%1.3f',v);
    else
        s = sprintf('%1.2f',v);
    end
else
    v = abs(v);
    if(v<0.001)
        s = '**';
    elseif(v<0.005)
        s = '*';
    elseif(v<0.01)
        s = '..';
    elseif(v<0.05)
        s = '.';
    else
        s = '';
    end
end
            
end