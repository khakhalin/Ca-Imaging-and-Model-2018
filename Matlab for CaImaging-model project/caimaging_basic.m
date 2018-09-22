function S = caimaging_basic(S,refCell)
% S = CAIMAGING_BASIC(S)
% S = CAIMAGING_BASIC(S,REFCELL)
%
% Finds spikes. It uses S.dataF and S.time to calculage S.dataS and S.timeS.
% Essentially it's my wrapper for fast_oopsi by Joshua Vogelstein.
% If REFCELL number is provided, it uses this cell to estimate noise parameters. You may use caimaging_browser()
% to find a good dead cell with proper noise level.

% Oct 26 2013: Created.
% Oct 31 2013: Updated with some more stuff.
% Dec 12 2013: Further development (towards something more meaningful).
% Dec 17 2013: fastoopsi integrated. Now this function can be considered a loading / pre-processing
%   function for further 'clever processing'.
% Dec 18 2013: It turned out that the way I exported stuff from NIS elements was incorrect (it supported
%   only 74 ROIs max), so now I need to re-design the reader to support another layout of input files.
% Jan 29 2014: + ActiveX way of learning number of sheets.
% Jan 31 2014: + Zero cell for noise calculations.
%               Reading removed to a separate function CAIMAGING_READ; only processing is left here.
% Mar 03 2014: + REFCELL parameter that allows to fix A and SIGMA options for fast_oopsi.


opt1 = [];                                              % Options structure placeholders.
opt2 = [];

fprintf('Thinking: ');
for(q=1:length(S))
    fprintf('.');
    
    dataF = smooth(S(q).time,S(q).dataF);               % Interpolate over all possible gaps (in a flat way, not even linear)
    [n,m] = size(dataF);
    dt = (max(S(q).time)-min(S(q).time))/n;
    S(q).timeS = (1:n)'*dt;   % Fix gaps in time as well
    
    opt1.dt = dt;                                       % Framerate (required)

    % zeroLevel = min(dataF);         % F0 for all subsequent (F-F0)/F0 calculations
    % zeroLevel = mean(dataF(1:10,:));
    % dataF = bsxfun(@plus, dataF,-zeroLevel);
    % dataF = bsxfun(@times,dataF,1./zeroLevel);

    %if(isfield(S,'noiseSigma'))
    if(nargin>1)                                        % Reference zero cell was provided
        if(~isfield(opt2,'sig'))
            fprintf('Estimating SIGMA and A from a reference cell %d\n',refCell);
            [~,P,~,~] = fast_oopsi(dataF(:,refCell)',opt1,opt2);
            opt2.sig = P.sig;
            opt2.a   = P.a;
            opt1.est_sig = 0;
            opt1.est_a = 0;
        end    
    end
    
    S(q).dataS = zeros(size(dataF));
    for(iCell=1:m)
        S(q).dataS(:,iCell) = fast_oopsi(dataF(:,iCell)',opt1,opt2)';
    end    
end
fprintf(' Done.\n');

% [cPca,sPca,eigenvalues] = princomp(amps);
% [~,i] = sort(sPca(:,1));
% figure; myplot(amps(i,:));
% ylabel('Cell number'); xlabel('Trial #');
% title('Response intensity (1st component)');
% 
% figure; plot(sPca(:,1),sPca(:,2),'bo');
% xlabel('Component 1'); ylabel('Component 2');
% title('Something like stimulus selectivity')

end



function b = myfilter(a,nWindow)

[n,m] = size(a);
a = [repmat(mean(a(1:nWindow,:)),nWindow,1) ; a ; repmat(mean(a(end-nWindow:end,:)),nWindow,1)];
[fb,fa] = butter(3,1/nWindow);
a = filtfilt(fb,fa,a);
b = a(nWindow+1:end-nWindow,:);

end


function b = smooth(t,a)
% Interpolate A over whatever gaps in T
jumps = [0; t(2:end)-t(1:end-1)];
% figure; plot(jumps)
step = mode(jumps);                   % Typical time step
[n,m] = size(a);
realN = round((max(t)-min(t))/step);
b = zeros(realN,m);
b(1,:) = a(1,:);
icounter = 1;

for(q=2:n)                              % Cycling through initial data file
    if(jumps(q)<1.5*step)               % Normal time step
        icounter = icounter+1;
        b(icounter,:) = a(q,:);
    else                                % A gap
        toAdd = round(jumps(q)/step);
        for(iAdd=1:toAdd)
            b(icounter+iAdd,:) = a(q,:) + (a(q+1,:)-a(q,:))/(toAdd+1)*iAdd;
        end
        icounter = icounter+toAdd;        
    end        
end

% for(q=2:realN)    
%     if(jumps(icounter)<1.5*step)   % Jumpt to this time tick was a normal single step
%         icounter = icounter+1;
%         b(q,:) = a(icounter,:);
%     else                           % A jump that needs to be smoothed away
%         jumps(icounter) = jumps(icounter)-step;
%         b(q,:) = a(icounter,:);
%     end    
% end 
end