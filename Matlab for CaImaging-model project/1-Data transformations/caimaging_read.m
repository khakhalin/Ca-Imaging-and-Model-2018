function res = caimaging_read2
% CAIMAGING_READ
%
% Reads Ca imaging data, and returns it as an output.
% Should explicitly (script-style) set the file name and the protocol key in the body of the function.

% Oct 26 2013:  Created.
% Oct 31 2013:  Updated with some more stuff.
% Dec 12 2013:  Further development (towards something more meaningful).
% Dec 17 2013:  fastoopsi integrated. Now this function can be considered a loading / pre-processing
%                   function for further 'clever processing'.
% Dec 18 2013:  It turned out that the way I exported stuff from NIS elements was incorrect (it supported
%                   only 74 ROIs max), so now I need to re-design the reader to support another layout of input files.
% Jan 29 2014:  + ActiveX way of learning number of sheets.
% Jan 31 2014:  + Zero cell for noise calculations.
%                   Now only reads data (as it is so slow). Next steps of signal processing are all done in CAIMAGING_BASIC function.
% Mar 11 2014:  + key assigned directly, and stored in S(1).key
% Apr 29 2014:  Minor updates.
% Mar 13 2015:  Superceded by version 2 (caimaging_read2)
% Mar 20 2017:  Now reads Excel files from one folder, and saves mat files to a different folder.


res = [];

%baseFolderIn = 'C:\_Data\___Ca imaging\_caimg s49\';
%baseFolderOut = 'C:\_Data\___Ca imaging\_caimg s49 mat\';

baseFolderIn = 'C:\_Data\___Ca imaging\_caimg s46\';
baseFolderOut = 'C:\_Data\___Ca imaging\_caimg s46 mat\';

iFolder = 0;
%%%% ---------- stage 46 set -------------
% iFolder = iFolder+1; folderName{iFolder} = '140718a';
% iFolder = iFolder+1; folderName{iFolder} = '140716b';
% iFolder = iFolder+1; folderName{iFolder} = '140716a';
% iFolder = iFolder+1; folderName{iFolder} = '140715';
% iFolder = iFolder+1; folderName{iFolder} = '140711';
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

%%%% ---------- stage 49 set -------------
%iFolder = iFolder+1; folderName{iFolder} = '140726';
%iFolder = iFolder+1; folderName{iFolder} = '140724';
%iFolder = iFolder+1; folderName{iFolder} = '140723';
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
% iFolder = iFolder+1; folderName{iFolder} = '140325';
% iFolder = iFolder+1; folderName{iFolder} = '140318';
% iFolder = iFolder+1; folderName{iFolder} = '140317';
% iFolder = iFolder+1; folderName{iFolder} = '140314';
% iFolder = iFolder+1; folderName{iFolder} = '140312';
% iFolder = iFolder+1; folderName{iFolder} = '140311';
% iFolder = iFolder+1; folderName{iFolder} = '140310'; % Strange

key = 'cfs';

Excel = actxserver('Excel.Application');        % Initialize excel server
Excel.DisplayAlerts = 0;

for(iFolder = 1:length(folderName));
    %%% -------------- Preprocessing -----------------
    if(folderName{iFolder}(end)~='\')
        folderName{iFolder} = [folderName{iFolder} '\'];
    end
    fullFolderName = [baseFolderIn folderName{iFolder}];
    fprintf('\nFolder: %s\n',fullFolderName);

    filesList = dir(fullFolderName);                           % Find all files in this directory
    fileName = [];                                                  % Flush
    counter = 0;
    if(length(filesList)>2)
        for(q=3:length(filesList))
            if(strcmp(filesList(q).name(end-3:end),'xlsx'))         % If Excel file
                counter = counter + 1;
                fileName{counter} = filesList(q).name;
                % fprintf(' File found: %s\n',fileName{counter});
            end
        end
    end
    nFiles = length(fileName);

    %%% -------------- Reading files -----------------
    for(iFile=1:nFiles)                             % For every file
        if(fileName{iFile}(1)=='~')                 % If temp file - skip the rest
            fprintf(' Skipping file (it''s not real): %s\n',fileName{iFile});  
            continue;                               
        end
        fprintf(' Opening file: %s\n',fileName{iFile});    
        ExcelWorkbook = Excel.workbooks.Open([fullFolderName fileName{iFile}]);       % This line takes pretty long to run, but not longer than trying to read an unexistent page
        nSheets = Excel.ActiveWorkbook.Sheets.Count;
        % nSheets = 12;   % Shortcut (unsafe if used blindly)

        for(iSheet=1:nSheets)            
            fprintf('  Reading sheet %d ... ',iSheet);

            ExcelWorkbook.Sheets.Item(iSheet).Activate();        
            Range = ExcelWorkbook.ActiveSheet.UsedRange;
            rawData = Range.Value;
            fprintf('Done. ');

            [data,text] = parse_data(rawData);

            [S(iSheet).time,S(iSheet).dataF,S(iSheet).xy] = deploy(data);   % Sort out this wonderful file layout
            S(iSheet).time = S(iSheet).time/1000;                           % From ms to seconds                
            S(iSheet).q = iSheet;   

            %matFileName = [folderName fileName{iFile}(1:end-5) '_' sprintf('%0d',iSheet) '.mat'];
            if(~exist([baseFolderOut folderName{iFolder}],'dir'))               % If the output folder doesn't exist, create it
                fprintf('folder %s does not exist. Creating.\n',folderName{iFolder});
                mkdir([baseFolderOut folderName{iFolder}]);                
            end
            matFileName = [baseFolderOut folderName{iFolder} fileName{iFile}(1:end-5) '.mat'];     % Update same file again and again
            fprintf('Writing mat-file: %s',matFileName');
            %save(matFileName,'data','text');
            save(matFileName,'S');
            fprintf(' .\n'); 
        end
        ExcelWorkbook.Close(false);
        fprintf(' File closed.\n');
    end
end
Excel.Quit;
Excel.delete;

end


function [t,d,xy] = deploy(data)
% T for time, D for data, XY for coordinates
% Columns:
% 1     2          3                4                5     6  7
% time, intensity, noiseEstimation, something silly, area, x, y

bordUp   = find(~isnan(data(2:end,1)) &  isnan(data(1:end-1,1)))+1; % When i-1 is NaN, while i is not (first row of each block)
bordDown = find( isnan(data(2:end,1)) & ~isnan(data(1:end-1,1)));   % When i+1 is NaN, while i is not (last row of each block)

bordUp   = [1; bordUp];
bordDown = [bordDown; length(data)];

nCells = length(bordUp);

t = data(bordUp(1):bordDown(1),1);
d = zeros(length(t),nCells);

for(q=1:nCells)
    d(:,q) = data(bordUp(q):bordDown(q),2);
    xy(q,1:2) = [data(bordUp(q),6:7)];
end

end



function [numericArray,textArray] = parse_data(data)
% PARSE_DATA parse data from raw cell array into a numeric array and a text
% cell array.
% [numericArray,textArray] = parse_data(data)
% Input:
%        data: cell array containing data from spreadsheet
% Return:
%        numericArray: double array containing numbers from spreadsheet
%        textArray: cell string array containing text from spreadsheet
%==========================================================================
%
% Taken from the Internet, not written by AKh (Mar 20 2015)

% ensure data is in cell array
if ischar(data)
    data = cellstr(data);
elseif isnumeric(data) || islogical(data)
    data = num2cell(data);
end

% Check if raw data is empty
if isempty(data)
    % Abort when all data cells are empty.
    textArray = {};
    numericArray = [];
    return
else
    % Trim empty leading and trailing rows
    % find empty cells
    emptycells = cellfun('isempty',data);
    nrows = size(emptycells,1);
    firstrow = 1;
    % find last of leading empty rows
    while (firstrow<=nrows && all(emptycells(firstrow,:)))
         firstrow = firstrow+1;
    end
    % remove leading empty rows
    data = data(firstrow:end,:);
    
    % find start of trailing empty rows
    nrows = size(emptycells,1);
    lastrow = nrows;
    while (lastrow>0 && all(emptycells(lastrow,:)))
        lastrow = lastrow-1;
    end
    % remove trailing empty rows
    data = data(1:lastrow,:);
    
    % find start of trailing NaN rows
    warning('off', 'MATLAB:nonIntegerTruncatedInConversionToChar');
    while (lastrow>0 && ~(any(cellfun('islogical', data(lastrow,:)))) && ...
                        all(isnan([data{lastrow,:}])))
        lastrow = lastrow-1;
    end
    warning('on', 'MATLAB:nonIntegerTruncatedInConversionToChar');
    % remove trailing NaN rows    
    data=data(1:lastrow,:);
    
    [n,m] = size(data);
    textArray = cell(size(data));
    textArray(:) = {''};
end

vIsNaN = false(n,m);

% find non-numeric entries in data cell array
vIsText = cellfun('isclass',data,'char');
vIsNaN = cellfun('isempty',data)|strcmpi(data,'nan')|cellfun('isclass',data,'char');

% place text cells in text array
if any(vIsText(:))
    textArray(vIsText) = data(vIsText);
else
    textArray = {};
end
% Excel returns COM errors when it has a #N/A field.
textArray = strrep(textArray,'ActiveX VT_ERROR: ','#N/A');

% place NaN in empty numeric cells
if any(vIsNaN(:))
    data(vIsNaN)={NaN};
end

% extract numeric data
data = reshape(data,n,m);
rows = size(data,1);
m = cell(rows,1);
% Concatenate each row first
for n=1:rows
    m{n} = cat(2,data{n,:});
end
% Now concatenate the single column of cells into a matrix
numericArray = cat(1,m{:});

    
% trim all-NaN leading rows and columns from numeric array
% trim all-empty trailing rows and columns from text arrays
[numericArray,textArray]=trim_arrays(numericArray,textArray);

% ensure numericArray is 0x0 empty.
if isempty(numericArray)
    numericArray = [];
end
end


function [numericArray,textArray] = trim_arrays(numericArray,textArray)
% trim leading rows or cols
% if the string result has dimensions corresponding to a column or row of
% zeros in the matrix result, trim the zeros.
%
% Source: http://www.mathworks.com/matlabcentral/answers/uploaded_files/1748/MODxlsread.m
if ~isempty(numericArray) && ~isempty(textArray)
    [mn, nn] = size(numericArray);
    [ms, ns] = size(textArray);

    if ms == mn
        % trim leading column(textArray) from numeric data
        firstcolm = 1;
        while (firstcolm<=nn && all(isnan(numericArray(:,firstcolm))))
            firstcolm = firstcolm+1;
        end
        numericArray=numericArray(:,firstcolm:end);
    end

    if ns == nn
        % trim leading NaN row(s) from numeric data
        firstrow = 1;
        while (firstrow<=mn && all(isnan(numericArray(firstrow,:))))
            firstrow = firstrow+1;
        end
        numericArray=numericArray(firstrow:end,:);
        
        % trim leading empty rows(s) from text data
        firstrow = 1;
        while (firstrow<=ms && all(cellfun('isempty',textArray(firstrow,:))))
            firstrow = firstrow+1;
        end
        textArray=textArray(firstrow:end,:);
    end
    
    % trim all-empty-string trailing rows from text array
	lastrow = size(textArray,1);
    while (lastrow>0 && all(cellfun('isempty',textArray(lastrow,:))))
        lastrow = lastrow-1;
    end
	textArray=textArray(1:lastrow,:);
    
    % trim all-empty-string trailing columns from text array
	lastcolm = size(textArray,2);
    while (lastcolm>0 && all(cellfun('isempty',textArray(:,lastcolm))))
        lastcolm = lastcolm-1;
    end
	textArray=textArray(:,1:lastcolm);

    % trim all-NaN trailing rows from numeric array
	lastrow = size(numericArray,1);
    while (lastrow>0 && all(isnan(numericArray(lastrow,:))))
        lastrow=lastrow-1;
    end
	numericArray=numericArray(1:lastrow,:);
    
    % trim all-NaN trailing columns from numeric array
	lastcolm = size(numericArray,2);
    while (lastcolm>0 && all(isnan(numericArray(:,lastcolm))))
        lastcolm=lastcolm-1;
    end
	numericArray=numericArray(:,1:lastcolm);
end
end