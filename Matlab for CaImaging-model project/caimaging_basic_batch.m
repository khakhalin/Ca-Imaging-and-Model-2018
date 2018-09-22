function caimaging_basic_batch()
% caimaging_basic_batch()
%
% Script-like funciton that relies on caimaging_basic and caimaging_find_good_cell,
% but works automatically and goes through all mat-files in the folder.
% Mat-files are supposedly created by caimaging_read2().
% Folder (or rather a sequence of them if necessary) is hard-coded in the function.

% Mar 20 2015: Created
% Mar 20 2017: Updated to work with a different folder.

res = [];
iFolder = 0;

%%% --- stage 49 set
% folderBaseIn = 'C:\_Data\___Ca imaging\_caimg s49 mat\';
% iFolder = iFolder+1; folderName{iFolder} = '140726';
% iFolder = iFolder+1; folderName{iFolder} = '140724';
% iFolder = iFolder+1; folderName{iFolder} = '140723';
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
% iFolder = iFolder+1; folderName{iFolder} = '140325'; % Some cell duplications
% iFolder = iFolder+1; folderName{iFolder} = '140318';
% iFolder = iFolder+1; folderName{iFolder} = '140317';
% iFolder = iFolder+1; folderName{iFolder} = '140314';
% iFolder = iFolder+1; folderName{iFolder} = '140312'; % Some cell duplications
% iFolder = iFolder+1; folderName{iFolder} = '140311'; % Some cell duplications

%%% --- stage 46 set
folderBaseIn = 'C:\_Data\___Ca imaging\_caimg s46 mat\';
iFolder = iFolder+1; folderName{iFolder} = '140718a';
iFolder = iFolder+1; folderName{iFolder} = '140716b';
iFolder = iFolder+1; folderName{iFolder} = '140716a';
iFolder = iFolder+1; folderName{iFolder} = '140715';
iFolder = iFolder+1; folderName{iFolder} = '140711';
iFolder = iFolder+1; folderName{iFolder} = '140708b';
iFolder = iFolder+1; folderName{iFolder} = '140708a';
iFolder = iFolder+1; folderName{iFolder} = '140705a';
iFolder = iFolder+1; folderName{iFolder} = '140704b';
iFolder = iFolder+1; folderName{iFolder} = '140704a';
iFolder = iFolder+1; folderName{iFolder} = '140627';
iFolder = iFolder+1; folderName{iFolder} = '140626';
iFolder = iFolder+1; folderName{iFolder} = '140620';
iFolder = iFolder+1; folderName{iFolder} = '140619b';
iFolder = iFolder+1; folderName{iFolder} = '140619a';
iFolder = iFolder+1; folderName{iFolder} = '140613';
iFolder = iFolder+1; folderName{iFolder} = '140611';
iFolder = iFolder+1; folderName{iFolder} = '140610';
iFolder = iFolder+1; folderName{iFolder} = '140530b';
iFolder = iFolder+1; folderName{iFolder} = '140530a';
iFolder = iFolder+1; folderName{iFolder} = '140529';
iFolder = iFolder+1; folderName{iFolder} = '140528';
iFolder = iFolder+1; folderName{iFolder} = '140516';
iFolder = iFolder+1; folderName{iFolder} = '140502';


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
    
    for(iFile=1:nFiles)
        S = [];                                                             % Flush
        matFileName = [folderBaseIn folderName{iFolder} fileName{iFile}];   % Update same file again and again
        fprintf('Reading    mat-file: %s ...',matFileName');
        load(matFileName,'S');
        fprintf(' Done.\n');
        S = checkAndFix(S);                                                 % Check for internal consistency and try to fix if necessary
        refCell = caimaging_find_good_cell(S,1);
        drawnow;        
        fprintf('Processing mat-file: %s\n',matFileName');
        S = caimaging_basic(S,refCell);
        fprintf(' Done.\n');
        %fprintf('Writing    mat-file: %s ...',matFileName');
        save(matFileName,'S');
        %fprintf(' Done.\n');        
    end
end

end


function S = checkAndFix(S)

nSweeps = length(S);
[nTime,nCells] = size(S(1).dataF);
for(is=2:nSweeps)
    nTime =  min(nTime, size(S(is).dataF,1));
    nCells = min(nCells,size(S(is).dataF,2));    
end

for(is=1:nSweeps)
    %fprintf('%3d ',size(S(is).dataF,2));
    if(size(S(is).dataF,2)>nCells)              % Unexpected number of cells        
        warning(sprintf('Sweep %d has a wrong number of cells (%d instead of %d). Using the first %d only',...
            is , size(S(is).dataF(1:nTime,:),2) , nCells, nCells));        
        S(is).dataF = S(is).dataF(:,1:nCells);
    end    
end
%fprintf('\n');

end