function caimaging_extra_figures3
% caimaging_extra_figures3
%
% This one is technical and makes the image of visual stimuli.

n = 120;        % Image resolution
ng = 5;         % Grid size for permutations
ns = n/ng;      % I trust that it's an integer; everything will break if it's not
[x,y] = meshgrid(1:n,1:n);
figure('Color','white');
for(iType = 1:3)
    p = randperm(ng^2);
    for(t=1:5)
        pic = zeros(n);
        d = n/2/4*(t-1);
        c = n/2;    % Center
        if(iType>1)
            pic((x-c).^2 + (y-c).^2 <= d^2) = 1;
        else
            pic(:) = 1*(t>2);
        end
        if(iType==3)
            temp = zeros(size(pic));
            for(i=1:ng)
                for(j=1:ng)
                    k = (i-1)*ng + j;
                    pi = floor((p(k)-1)/ng)+1;
                    pj = p(k)-(pi-1)*ng;
                    temp((1:ns)+ns*(i-1),(1:ns)+ns*(j-1)) = pic((1:ns)+ns*(pi-1),(1:ns)+ns*(pj-1));
                end
            end
            pic = temp;
        end
        pic(:,[1 2 end-1 end]) = 0; pic([1 2 end-1 end],:) = 0;  % Frame
        subplot(3,5,(iType-1)*5+t);        
        myplot(pic);
        set(gca,'Visible','off','PlotBoxAspectRatio',[1 1 1]);
        colormap('gray'); caxis([0 1]);
        drawnow();
    end
end

end