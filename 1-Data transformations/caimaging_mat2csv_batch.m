function caimaging_mat2csv_batch()
% caimaging_mat2cvs_batch()
%
% Script-like funciton that reads ca imaging mat files and saves them again as
% csv (comma-separated values) files. Assumes that mat files were created by
% a sequence of caimaging_read2, with post-processing by caimaging_basic_batch.

% Mar 20 2017: Created. 

res = [];

stage = 46;

if(stage==49)
    folderBaseIn = 'C:\_Data\___Ca imaging\_caimg s49 mat\';
    folderBaseOut = 'C:\_Data\___Ca imaging\_caimg s49 csv\';

    iFolder = 0;
    iFolder = iFolder+1; folderName{iFolder} = '140726';
    % % iFolder = iFolder+1; folderName{iFolder} = '140724';
    % % iFolder = iFolder+1; folderName{iFolder} = '140723';
    % iFolder = iFolder+1; folderName{iFolder} = '140722';
    % iFolder = iFolder+1; folderName{iFolder} = '140718b';
    % iFolder = iFolder+1; folderName{iFolder} = '140714';
    % iFolder = iFolder+1; folderName{iFolder} = '140710';
    % iFolder = iFolder+1; folderName{iFolder} = '140709';
    % iFolder = iFolder+1; folderName{iFolder} = '140707';
    % iFolder = iFolder+1; folderName{iFolder} = '140705b';
    % iFolder = iFolder+1; folderName{iFolder} = '140703';
    % iFolder = iFolder+1; folderName{iFolder} = '140612';
    % iFolder = iFolder+1; folderName{iFolder} = '140522';
    % iFolder = iFolder+1; folderName{iFolder} = '140521';
    % iFolder = iFolder+1; folderName{iFolder} = '140505';
    % iFolder = iFolder+1; folderName{iFolder} = '140408';
    % iFolder = iFolder+1; folderName{iFolder} = '140328';
    % iFolder = iFolder+1; folderName{iFolder} = '140326';
    % % iFolder = iFolder+1; folderName{iFolder} = '140325'; %%% Some problems with this experiment - it generates an error in cfs3
    % iFolder = iFolder+1; folderName{iFolder} = '140318';
    % iFolder = iFolder+1; folderName{iFolder} = '140317';
    % iFolder = iFolder+1; folderName{iFolder} = '140314';
    % % iFolder = iFolder+1; folderName{iFolder} = '140312'; %%% Some problems with this experiment - it generates an error in cfs1
    % % iFolder = iFolder+1; folderName{iFolder} = '140311'; %%% Some problems with this experiment - it generates an error in cfs1
else
    folderBaseIn = 'C:\_Data\___Ca imaging\_caimg s46 mat\';
    folderBaseOut = 'C:\_Data\___Ca imaging\_caimg s46 csv\';
    %%% Good recordings only; weak are super-commented out

    iFolder = 0;
    %%% iFolder = iFolder+1; folderName{iFolder} = '140502';
    iFolder = iFolder+1; folderName{iFolder} = '140516';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140528';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140529';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140530a';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140530b';
    iFolder = iFolder+1; folderName{iFolder} = '140610';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140611';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140613';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140619a';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140619b';
    %%% iFolder = iFolder+1; folderName{iFolder} = '140620';
    iFolder = iFolder+1; folderName{iFolder} = '140626';
    iFolder = iFolder+1; folderName{iFolder} = '140627';
    iFolder = iFolder+1; folderName{iFolder} = '140704a';
    iFolder = iFolder+1; folderName{iFolder} = '140704b';
    iFolder = iFolder+1; folderName{iFolder} = '140705a';
    iFolder = iFolder+1; folderName{iFolder} = '140708a';
    iFolder = iFolder+1; folderName{iFolder} = '140708b';
    iFolder = iFolder+1; folderName{iFolder} = '140711';
    iFolder = iFolder+1; folderName{iFolder} = '140715';
    iFolder = iFolder+1; folderName{iFolder} = '140716a';
    iFolder = iFolder+1; folderName{iFolder} = '140716b';
    iFolder = iFolder+1; folderName{iFolder} = '140718a';    
end

key = 'cfs';

for(iFolder = 1:length(folderName));
    %%% -------------- Preprocessing -----------------
    if(folderName{iFolder}(end)~='\')
        folderName{iFolder} = [folderName{iFolder} '\'];
    end
    fullFolderName = [folderBaseIn folderName{iFolder}];
    fprintf('\nFolder: %s\n',fullFolderName);

    filesList = dir(fullFolderName);                                        % Find all files in this directory
    fileName = [];                                                          % Flush memory
    counter = 0;
    if(length(filesList)>2)
        for(q=3:length(filesList))
            if(strcmp(filesList(q).name(end-2:end),'mat'))                  % If matlab data file
                counter = counter + 1;
                fileName{counter} = filesList(q).name;
                % fprintf(' File found: %s\n',fileName{counter});
            end
        end
    end
    nFiles = length(fileName);
    
    Fbag = [];          % Here fluorescense data will be combined   
    Sbag = [];          % Here spiking esitmations will be combined
    sweepCounter = 1;
    for(iFile=1:nFiles)                 
        S = [];                                                             % Flush, just in case
        matFileName = [folderBaseIn folderName{iFolder} fileName{iFile}];
        fprintf('Reading    mat-file: %s\n',matFileName');
        load(matFileName,'S');
        for(iSweep=1:length(S))
            Fbag = [Fbag; ones(size(S(iSweep).time))*sweepCounter S(iSweep).time S(iSweep).dataF];
            Sbag = [Sbag; ones(size(S(iSweep).timeS))*sweepCounter S(iSweep).timeS S(iSweep).dataS];
            sweepCounter = sweepCounter+1;
        end        
    end
    
    folderOut = [folderBaseOut folderName{iFolder}];
    if(~exist(folderOut,'dir'))               % If the output folder doesn't exist, create it
        fprintf('folder %s does not exist. Creating.\n',folderName{iFolder});
        mkdir(folderOut);                
    end
    
    fprintf('Writing experiment %s: ',folderName{iFolder});
    
    fid = fopen([folderOut 'readme.txt'],'w');
    fprintf(fid,'Recording id: %s\r\n',folderName{iFolder});
    fprintf(fid,'N cells: %d\r\n',size(S(1).xy,1));
    fprintf(fid,'F file for fluoresence (raw data); S file for spike estimations (deconvolved).\r\n');
    fprintf(fid,'In both files 1st column is sweep number, 2nd column is time stamp, remaining columns are data.\r\n');
    fprintf(fid,'Stimuli always cycle through three stimulus types: looming stimulus, full-field flash, scrambled stimulus.\r\n');
    fclose(fid);
    
    temp = Fbag(:,3:end);
    Fbag(:,3:end) = round(temp);    % It seems that most neurons are whole numbers, but some have 1 sign after point, and I don't think it's needed.
    
    temp = Sbag(:,3:end);
    temp(temp<0.0001) = 0;          % To get rid of 1e-11 that are probably useless    
    Sbag(:,3:end) = temp;
    
    figure;                 % Debugging figure
    subplot(1,2,1); hold on;
    subplot(1,2,2); hold on;
    for(iSweep = 1:max(Sbag(:,1)))
        subplot(1,2,1); 
        g = find(Fbag(:,1)==iSweep);    % Logical subset. Find is needed for correct plotting ofset in the next row.
        plot(Fbag(g,2) + 4.1*mod(iSweep-1,3) , mean(Fbag(g,3:end),2)-mean(Fbag(g(1),3:end))+200*floor((iSweep-1)/3));
        subplot(1,2,2); 
        g = Sbag(:,1)==iSweep;          % Logical subset
        plot(Sbag(g,2) + 4.1*mod(iSweep-1,3) , Sbag(g,3:end)+floor((iSweep-1)/3));        
    end
    title(folderName{iFolder}(1:end-1)); 
    drawnow();
    
    dlmwrite([folderOut 'xy.txt'],S(1).xy,'delimiter',',','newline','pc');
    dlmwrite([folderOut sprintf('%sF.txt',folderName{iFolder}(1:end-1))],Fbag,'delimiter',',','newline','pc');
    dlmwrite([folderOut sprintf('%sS.txt',folderName{iFolder}(1:end-1))],Sbag,'delimiter',',','newline','pc');
    fprintf('Finished.\n');
end

end