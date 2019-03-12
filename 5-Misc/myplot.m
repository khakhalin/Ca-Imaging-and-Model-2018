function h = myplot(varargin);
% h = myplot(Data)
% h = myplot(AX, Data)
% h = myplot(AX, Data, 'colormap') 
% h = myplot(X, Y, Data, 'colormap') 
% h = myplot(AX, X, Y, Data, 'colormap') 
%
% Draws a Heatmap, returns the surf object handle.
% If axes handle AX is given, outputs to it.
% 'hot' colormap is used by default

useMap = 'hot'; % Default value

[cax,ars,nargs] = axescheck(varargin{:});

if(isempty(cax)) % Axes were not provided
    cax = gca;
    arguments = varargin(1:end);
else
    arguments = varargin(2:end);
end

switch length(arguments)
    case 1
        Data = arguments{1};
    case 2
        Data = arguments{1};
        useMap = varargin{2};        
    case 3
        x = arguments{1};
        y = arguments{2};
        Data = arguments{3};        
    case 4
        x = arguments{1};
        y = arguments{2};
        Data = arguments{3};
        useMap = varargin{4};        
end        
h = surf((1:size(Data,2)+1)-.5, (1:size(Data,1)+1)-.5, zeros(size(Data,1)+1, size(Data,2)+1), Data);        

set(gca,'XLim',[0.5 size(Data,2)+0.5]);
set(gca,'YLim',[0.5 size(Data,1)+0.5]);
view(cax,2);
if(length(arguments)>=3)
    set(gca,'XTick',1:length(x),'XTickLabel',x);
    set(gca,'YTick',1:length(y),'YTickLabel',y);
end
switch useMap
    case 'abone'
        colormap(cax,'bone');
        useMap = get(gcf,'ColorMap');
        useMap = useMap(end:-1:1,:);
    case 'ahot'
        colormap('hot');
        useMap = 1-get(gcf,'ColorMap');
    otherwise
        % Nothing
end
colormap(cax,useMap);
shading(cax,'flat');

end
