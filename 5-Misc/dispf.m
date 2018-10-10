    function dispf(data,format)
% DISPF(DATA)
% DISPF(DATA,FORMAT)
%
% Works kind of in between DISP and FPRINTF: prints the data rectangularly, but can follow formats.
% 
% Example: dispf(a,'$5f');

% May 06 2013: Updated.
% May 07 2013: Partial cell support.
% May 15 2014: + simple strings are just returned as strings.
% Sep 13 2014: Better support of strings.


[n,m] = size(data);
if(nargin<2)
    if(~ischar(data))
        format = '%5f';
    else
        fprintf([data '\n']);
        return;
    end
end

if(~iscell(data))
    if(m>1)
        for(q=1:n)
            fprintf([format '\t'],data(q,1:end-1));
            fprintf([format '\n'],data(q,end));
        end
    else
        fprintf([format '\n'],data);
    end
else % Cell
    if(ischar(data{1}))
        if(nargin<2)
            format = '%s';
        end
        if(m>1)
            for(q=1:n)
                fprintf([format '\t'],data{q,1:end-1});
                fprintf([format '\n'],data{q,end});
            end
        else
            for(q=1:n)
                fprintf([format '\n'],data{q});
            end
        end
    else % Not char
        error('Error: DISPF does not yet know how to display numerical cell arrays.');
    end
end

end