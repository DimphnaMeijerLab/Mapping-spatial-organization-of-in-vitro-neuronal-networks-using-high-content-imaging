% Initialize and setup
delete(gcp('nocreate'));
clc
clearvars
close all

%% Add Module Paths
workfolder = pwd;

addpath(...
    fullfile(workfolder),...
    fullfile(workfolder, "Distance Analysis"),...
    fullfile(workfolder, 'Nuclei Analysis'),...
    fullfile(workfolder, 'Plotting'));

%% Data Selection

% Input path
disp('Select the folder that contains your Data')
datapath = uigetdir('E:\', 'Select the folder that contains your Data');

% List Datasets in Data folder
tempdatasets = GetSubDirsFirstLevelOnly(datapath);

% Select desired Folders/Datasets in Data folder
include = listdlg('ListString', tempdatasets, 'PromptString', ...
    'Please select the folder(s) that you want to analyse', 'ListSize', [300, 300]);

% Select Datasets
datasets = [];
for i = 1:length(tempdatasets)
    dataset = string(tempdatasets(i));
    if ismember(dataset, tempdatasets(include))
        datasets = [datasets; dataset]; %#ok<AGROW>
    end
end

% Output path
disp("Select the folder in which you'd like to save your Results")
savepath = uigetdir('E:\', "Select the folder in which you'd like to save your Results");

%% Enter Worker Count
workercount = 16;

% Initialize cell arrays to collect data
NtoNCell = cell(length(datasets), 1);
NtoCCell = cell(length(datasets), 1);
NtoNCCell = cell(length(datasets), 1);
NtoNLCell = cell(length(datasets), 1);

%----Parfor loop over Datasets----
tic
parpool('local', workercount);
cellResultData = cell(length(datasets), 3);

parfor i = 1:length(datasets) % parfor

    % Local cell arrays for this dataset
    localNtoN = [];
    localNtoC = [];
    localNtoNC = [];
    localNtoNL = [];

    dataset = string(datasets(i));
    disp("Processing Dataset: " + dataset);

    % List Wells, make string from cell array and transpose for data struct
    wells = string(GetSubDirsFirstLevelOnly(fullfile(datapath, dataset)))';

    cellResultWells = cell(length(wells), 2);

    %----Loop over Wells----
    for j = 1:length(wells)

        well = wells(j);
        disp('Processing well: ' + well + ' of ' + dataset);

        % Load neuron centroids from Excel
        % Assume the Excel file is named as {well}_neuron_centroids.xlsx
        excelFilePath = fullfile(datapath, dataset, well, well + "_neuron_centroids.xlsx");
        if isfile(excelFilePath)
            neuronCentroids = readmatrix(excelFilePath, 'Sheet', 'Sheet1');
        else
            disp(['Excel file not found for well ' well '. Skipping...']);
            continue;
        end

        % Ensure the centroids are in [x, y] format
        if size(neuronCentroids, 2) ~= 2
            disp(['Invalid centroid data in ' excelFilePath '. Skipping...']);
            continue;
        end

        % Call Distance Analysis
        disp("Calculating distances...");
        [NtoN, NtoC, NtoNC, NtoNL] = Distance_analysis_centroids_v3(neuronCentroids);
        disp("Distance analysis complete");

        % Save well results
        wellResults = struct(...
            "NeuronCentroids", neuronCentroids, ...
            'NtoN', NtoN, ...
            'NtoC', NtoC, ...
            'NtoNC', NtoNC, ...
            'NtoNL', NtoNL ...
        );

        cellResultWells(j, :) = {well, wellResults};

        % Save the well results
        wellSavePath = fullfile(savepath, dataset, well);
        if ~isfolder(wellSavePath)
            mkdir(wellSavePath);
        end
        parsave(fullfile(wellSavePath, well + ".mat"), wellResults);
        
        % Collect local distances row-wise
        localNtoN = [localNtoN; NtoN(:)'];
        localNtoC = [localNtoC; NtoC(:)'];
        localNtoNC = [localNtoNC; NtoNC(:)'];
        localNtoNL = [localNtoNL; NtoNL(:)'];
    end

    % Store local distances in the cell arrays
    NtoNCell{i} = localNtoN;
    NtoCCell{i} = localNtoC;
    NtoNCCell{i} = localNtoNC;
    NtoNLCell{i} = localNtoNL;

    % Save Dataset Results
    datasetResults = struct("well", num2cell(wells), "NeuronCentroids", {cellResultWells(:, 2)});
    cellResultData(i, :) = {dataset, datasetResults, cellResultWells};

    datasetSavePath = fullfile(savepath, dataset);
    if ~isfolder(datasetSavePath)
        mkdir(datasetSavePath);
    end
    parsave(fullfile(datasetSavePath, dataset + ".mat"), datasetResults);

    disp("Finished processing " + dataset)
end

toc
delete(gcp('nocreate'));
disp("Finished all datasets");
disp("Script Complete!");

% Combine all distances into a single struct, row-wise
allDistancesStruct = struct(...
    'NtoN', cat(1, NtoNCell{:}), ...
    'NtoC', cat(1, NtoCCell{:}), ...
    'NtoNC', cat(1, NtoNCCell{:}), ...
    'NtoNL', cat(1, NtoNLCell{:}) ...
);

% Save the struct with all distances
save(fullfile(savepath, 'distances.mat'), 'allDistancesStruct');

% parsave function to save variables inside parfor loop
function parsave(fname, var)
    save(fname, 'var', '-v7.3');
end