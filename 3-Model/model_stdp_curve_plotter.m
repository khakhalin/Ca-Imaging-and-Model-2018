function varargout = model_stdp_curve_plotter(type,oneFlag)
% model_stdp_curve_plotter
%
% Similar to model_stdp_tester(), in the sense that it reads files created by model_stdp().
% But this one does not run any simulations, and only builds the curves that are already saved
% in those model-output files. As of Apr 2019, these curves are not used in any figures in the paper.


%%% ------------------------------------ Constants ------------------------------------
inFolder = 'C:\Users\Arseny\Documents\3_Modeling\2018 Model outputs\';  % A folder to look for input files
outFile = ['C:\Users\Arseny\Documents\3_Modeling\modelAnalysis' datestr(now,'yymmdd') '.csv'];     % A file to put all output values into
if(nargin<1)
    type = 'all';                                   % What type of training files to collect
end
if(nargin<2); oneFlag = 0; end;                     % If not given, assume that we want full analysis, not just one plot
M = [];                                             % In this structure we will keep a table of all data. See functions PUSH, REMEMBER, and the saving block below

flagFigDeg = 0;                                     % Whether we want a figure with degree histograms or not

%%% -------- Read the data
fprintf('Started. Looking for file type: %s\n',type);
fileList = dir(inFolder);                       % List of files
fileList = {fileList(3:end).name};              % Skip '.' and '..' that start the list (at least on Windows)
if(oneFlag); fileList = fileList(1); end;       % If we need only one plot, let there be one plot

%%% -------- Initialize the cumulative figure object
f.hF = figure('color','white');
for(i=1:6)
    f.hp(i) = subplot(2,3,i); set(gca,'FontSize',8); hold on;
end
title(f.hp(1),'Efficiency');                %set(f.hp(1),'YLim',[0 0.1]);
title(f.hp(2),'Modularity');                set(f.hp(2),'YLim',[0 1]);
title(f.hp(3),'cluster');                                       set(f.hp(3),'YLim',[0 0.03]);
title(f.hp(4),'Flow');                                         %set(f.hp(4),'YLim',[0 10]);
title(f.hp(5),'antiFlow');                                     %set(f.hp(5),'YLim',[0 10]);

set(f.hF,'UserData',f);


%%% ------------ Main part
bagResults = [];                                % To keep all OUT values
bagMeasures = [];                               % To store network measures (those stored in data files), for the summary

for(iFile=1:length(fileList))
    try
        temp = load([inFolder fileList{iFile}]);
    catch
        fprintf('Unreadable file: %s\n',fileList{iFile});
        continue;
    end    
    U = temp.U;                                     % Fetch the data    
    if(strcmp(type,U.type) || strcmp(type,'all')) 	% Either correct experiment type or a wildcard
        % Explanation of U.measurements (1st column is time)
        % dataLabels = {'Efficiency','InInAssort','Modular','Cluster','Flow','revFlow'};
        %                2            3               4       5         6         7     
        bagMeasures = cat(3,bagMeasures,U.measurements);    % Concat in 3D: now time runs down, var type right, and experiment# deep
        plot(f.hp(1),U.measurements(:,1),U.measurements(:,2),'-','Color',[.7 1 .6]); % Efficiency
        plot(f.hp(2),U.measurements(:,1),U.measurements(:,4),'-','Color',[.7 1 .6]); % Modularity
        plot(f.hp(3),U.measurements(:,1),U.measurements(:,5),'-','Color',[.7 1 .6]); % Cluster
        plot(f.hp(4),U.measurements(:,1),U.measurements(:,6),'-','Color',[.7 1 .6]); % rank        
        plot(f.hp(5),U.measurements(:,1),U.measurements(:,7),'-','Color',[.7 1 .6]); % rev rank
    end
end

plot(f.hp(1),bagMeasures(:,1,1),mean(bagMeasures(:,2,:),3),'k-');   % Eff
plot(f.hp(2),bagMeasures(:,1,1),mean(bagMeasures(:,4,:),3),'k-');   % Mod
plot(f.hp(3),bagMeasures(:,1,1),mean(bagMeasures(:,5,:),3),'k-');   % cluster
plot(f.hp(4),bagMeasures(:,1,1),mean(bagMeasures(:,6,:),3),'k-');   % rank
plot(f.hp(5),bagMeasures(:,1,1),mean(bagMeasures(:,7,:),3),'k-');   % rev rank

plot(f.hp(1),bagMeasures(:,1,1),mean(bagMeasures(:,2,:),3)-std(bagMeasures(:,2,:),[],3),'b-');   % Eff
plot(f.hp(2),bagMeasures(:,1,1),mean(bagMeasures(:,4,:),3)-std(bagMeasures(:,4,:),[],3),'b-');   % Mod
plot(f.hp(3),bagMeasures(:,1,1),mean(bagMeasures(:,5,:),3)-std(bagMeasures(:,5,:),[],3),'b-');   % cluster
plot(f.hp(4),bagMeasures(:,1,1),mean(bagMeasures(:,6,:),3)-std(bagMeasures(:,6,:),[],3),'b-');   % rank
plot(f.hp(5),bagMeasures(:,1,1),mean(bagMeasures(:,7,:),3)-std(bagMeasures(:,7,:),[],3),'b-');   % rev rank

plot(f.hp(1),bagMeasures(:,1,1),mean(bagMeasures(:,2,:),3)+std(bagMeasures(:,2,:),[],3),'b-');   % Eff
plot(f.hp(2),bagMeasures(:,1,1),mean(bagMeasures(:,4,:),3)+std(bagMeasures(:,4,:),[],3),'b-');   % Mod
plot(f.hp(3),bagMeasures(:,1,1),mean(bagMeasures(:,5,:),3)+std(bagMeasures(:,5,:),[],3),'b-');   % cluster
plot(f.hp(4),bagMeasures(:,1,1),mean(bagMeasures(:,6,:),3)+std(bagMeasures(:,6,:),[],3),'b-');   % rank
plot(f.hp(5),bagMeasures(:,1,1),mean(bagMeasures(:,7,:),3)+std(bagMeasures(:,7,:),[],3),'b-');   % rev rank

return;

end
