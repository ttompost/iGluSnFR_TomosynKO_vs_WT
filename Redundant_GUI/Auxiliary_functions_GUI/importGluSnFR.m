function gluT = importGluSnFR()

% Get the folder with the data
dataDir = uigetdir(pwd);
filePattern = fullfile(dataDir, '*.csv');
dataFiles = dir(filePattern);
gluT = table;
tempT = table;
% Import the data
hWait = waitbar(0, 'Loading data');
for file = 1:numel(dataFiles)
    waitbar(file/numel(dataFiles), hWait, sprintf('Loading data (%.2f%%)', file/numel(dataFiles)*100));
    filePath = fullfile(dataFiles(file).folder, dataFiles(file).name);
    tempData = readtable(filePath);
    tempData = table2array(tempData(:,2:end));
    % Calculate the FF0 (condiser the first 5 frames as baseline)
    [nFrames, nRoi] = size(tempData);
    baseInt = mean(tempData(1:5,:));
    deltaff0Ints = (tempData - repmat(baseInt, nFrames, 1)) ./ repmat(baseInt, nFrames, 1);
    % Create the final table data
    fileName = regexp(dataFiles(file).name, '_', 'split');
    weekID = fileName{2};
    recID = fileName{3};
    covID = fileName{4};
    cellID = fileName{5}(1:3);
    if file == 1
        gluT.Filename = {[weekID '_' recID '_' covID '_' cellID]};
        gluT.week = {weekID};
        gluT.Recording = {recID};
        gluT.Coverslip = {covID};
        gluT.CellID = {cellID};
        gluT.RawData = {tempData'};
        gluT.FF0Data = {deltaff0Ints'};
    else
        tempT.Filename = {[weekID '_' recID '_' covID '_' cellID]};
        tempT.week = {weekID};
        tempT.Recording = {recID};
        tempT.Coverslip = {covID};
        tempT.CellID = {cellID};
        tempT.RawData = {tempData'};
        tempT.FF0Data = {deltaff0Ints'};
        gluT = [gluT; tempT];
%         if strcmp(recID, 'rec1')
%             gluT = [gluT; tempT];
%         else
%             sameCell = contains(gluT.Coverslip, covID) & contains(gluT.CellID, cellID);
%             gluT{sameCell, 'RawData'} = {[cell2mat(gluT{sameCell, 'RawData'}), nan(nRoi,1), tempData']};
%             gluT{sameCell, 'FF0Data'} = {[cell2mat(gluT{sameCell, 'FF0Data'}), nan(nRoi,1), deltaff0Ints']};
%         end
    end
end
close(hWait);
gluT.week = categorical(gluT.week);
gluT.Recording = categorical(gluT.Recording);
gluT.MaxData = cellfun(@max, gluT.FF0Data, 'UniformOutput', false);