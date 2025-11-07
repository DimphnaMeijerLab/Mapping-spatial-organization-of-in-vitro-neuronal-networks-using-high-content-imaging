%% MIA - Modular Image Analysis
% 
% Main analysis script for whole-well image datasets.
% This script:
%   - Iterates over selected datasets and wells
%   - Runs nuclei analysis (neurons / glia segmentation and stats)
%   - Runs distance analysis between neuronal and glial centroids
%   - Saves per-well and per-dataset results
%
% Requirements:
%   - Functions:
%       GetSubDirsFirstLevelOnly.m
%       Nuclei_analysis.m
%       Distance_analysis.m         % must accept pixelSize & wellRadius
%       parsave.m
%   - Folder structure:
%       datapath/
%           dataset1/
%               wellA/
%                   wellA_th_dapi.tif
%                   wellA_th_AnkG.tif
%                   wellA_dapi.tif
%               wellB/
%                   ...
%
% Notes:
%   - Distances are computed using Distance_analysis and can be converted
%     to physical units inside that function using the provided pixelSize.
%   - wellRadius is passed for normalisation / quality control, as needed.

%% Housekeeping
delete(gcp('nocreate'));
clc;
clearvars;
close all;

%% Add Module Paths
workfolder = pwd;

addpath( ...
    fullfile(workfolder), ...
    fullfile(workfolder, "Distance Analysis"), ...
    fullfile(workfolder, "Nuclei Analysis"), ...
    fullfile(workfolder, "Plotting") ...
    );

%% Imaging / analysis parameters
% Define imaging parameters here for reproducibility.
% Adjust these values according to your microscope / assay.

pixelSize  = 0.401;   % [µm / pixel]   - linear pixel size
wellRadius = 3250;   % [µm]           - physical radius of the well

%% Data selection

% Input path
disp('Select the folder that contains your Data');
datapath = uigetdir('E:\', ...
    'Select the folder that contains your Data');

if datapath == 0
    error('No data folder selected. Aborting.');
end

% List datasets in data folder (first-level subdirectories)
tempdatasets = GetSubDirsFirstLevelOnly(datapath);

% Select desired datasets
include = listdlg( ...
    'ListString',      tempdatasets, ...
    'PromptString',    'Please select the folder(s) that you want to analyse', ...
    'ListSize',        [300, 300] ...
    );

if isempty(include)
    error('No datasets selected. Aborting.');
end

datasets = string(tempdatasets(include));

% Output path
disp('Select the folder in which you''d like to save your Results');
savepath = uigetdir('E:\', ...
    'Select the folder in which you''d like to save your Results');

if savepath == 0
    error('No save folder selected. Aborting.');
end

%% Worker count for parallel processing
workercount = 16;

%% ---- Parfor loop over Datasets ----
tic;
parpool('local', workercount);

cellResultData = cell(length(datasets), 4);

parfor i = 1:length(datasets)

    dataset = string(datasets(i));
    disp("Processing Dataset: " + dataset);

    % List wells (first-level subdirectories of the dataset)
    wells = string(GetSubDirsFirstLevelOnly(fullfile(datapath, dataset)))';

    % Preallocate variables

    % Label images (logical)
    LBW_Dapi    = false(2);
    LBW_AnkG    = false(2);
    LBW_Neurons = false(2);

    % Nuclei analysis
    Neurons          = zeros(length(wells), 1);
    Neuron_stats     = struct();
    Neuron_stats_all = cell(1, length(wells));

    Glia             = zeros(length(wells), 1);
    Glia_stats       = struct();
    Glia_stats_all   = cell(1, length(wells));

    Thresholds       = zeros(length(wells), 1);

    % Distance analysis accumulators
    NtoN    = zeros(1);
    GtoG    = zeros(1);
    NtoC    = zeros(1);
    GtoC    = zeros(1);
    NtoNC   = zeros(1);
    GtoGC   = zeros(1);
    NtoNL   = zeros(1);
    GtoGL   = zeros(1);
    NtoGL   = zeros(1);
    GtoNL   = zeros(1);

    AllNtoN  = zeros(1);
    AllGtoG  = zeros(1);
    AllNtoC  = zeros(1);
    AllGtoC  = zeros(1);
    AllNtoNC = zeros(1);
    AllGtoGC = zeros(1);
    AllNtoNL = zeros(1);
    AllGtoGL = zeros(1);
    AllNtoGL = zeros(1);
    AllGtoNL = zeros(1);

    % Storage for per-well results
    cellResultWells = cell(length(wells), 2);

    %% ---- Loop over Wells ----
    for j = 1:length(wells)

        well = wells(j);
        disp("  Processing Well: " + well + " of " + dataset);

        % Build basic file prefix (e.g. 'A01')
        wellChar = char(well);

        %% Load Required Image Files

        % Thresholded DAPI (nuclei)
        LBW_Dapi = logical(imread(fullfile( ...
            datapath, dataset, well, [wellChar '_th_dapi.tif'])));

        % Thresholded AnkyrinG (axon initial segments)
        LBW_AnkG = logical(imread(fullfile( ...
            datapath, dataset, well, [wellChar '_th_AnkG.tif'])));

        % Raw DAPI image
        Dapi = imread(fullfile( ...
            datapath, dataset, well, [wellChar '_dapi.tif']));

        %% Module 1 – Nuclei Analysis
        disp("    Analysing nuclei...");

        [ ...
            Glia_stats, Neuron_stats, ...
            LBW_Glia, LBW_Neurons, ...
            Thresholds(j, 1) ...
        ] = Nuclei_analysis(LBW_Dapi, LBW_AnkG, Dapi);

        % Count neurons and glia
        Neurons(j) = length(Neuron_stats.Area);
        Glia(j)    = length(Glia_stats.Area);

        disp("    Nuclei analysis complete.");
        disp("    Calculating distances...");

        %% Module 2 – Distance Analysis
        % Distance_analysis should internally use pixelSize to convert
        % pixel distances to physical units (e.g. µm), and may use
        % wellRadius for normalisation / QC.

        [ ...
            NtoN, GtoG, NtoC, ...
            GtoC, NtoNC, GtoGC, ...
            NtoNL, GtoGL, NtoGL, GtoNL ...
        ] = Distance_analysis( ...
                Neuron_stats.Centroid, ...
                Glia_stats.Centroid, ...
                pixelSize, ...
                wellRadius ...
            );

        % Accumulate distances over wells
        AllNtoN  = [AllNtoN,  NtoN];
        AllGtoG  = [AllGtoG,  GtoG];
        AllNtoC  = [AllNtoC,  NtoC];
        AllGtoC  = [AllGtoC,  GtoC];
        AllNtoNC = [AllNtoNC, NtoNC];
        AllGtoGC = [AllGtoGC, GtoGC];
        AllNtoNL = [AllNtoNL, NtoNL];
        AllGtoGL = [AllGtoGL, GtoGL];
        % If needed, uncomment:
        % AllNtoGL = [AllNtoGL, NtoGL];
        % AllGtoNL = [AllGtoNL, GtoNL];

        disp("    Distance analysis complete.");
        disp("  Well " + well + " of " + dataset + " finished.");

        %% Save per-well results into struct

        wellResults = struct( ...
            "Neurons",     LBW_Neurons, ...
            "Glia",        LBW_Glia, ...
            "Neuronstats", Neuron_stats, ...
            "Gliastats",   Glia_stats, ...
            "NtoN",        NtoN, ...
            "GtoG",        GtoG, ...
            "NtoC",        NtoC, ...
            "GtoC",        GtoC, ...
            "NtoNC",       NtoNC, ...
            "GtoGC",       GtoGC, ...
            "NtoNL",       NtoNL, ...
            "GtoGL",       GtoGL, ...
            "NtoGL",       NtoGL, ...
            "GtoNL",       GtoNL ...
            );

        cellResultWells(j, :) = {well, wellResults};

        %% Save per-well results to disk

        wellSaveDir = fullfile(savepath, dataset, well);
        if ~isfolder(wellSaveDir)
            mkdir(wellSaveDir);
        end

        parsave(fullfile(wellSaveDir, [wellChar '.mat']), wellResults);

    end % wells loop

    %% Save Dataset Results (summary)

    datasetResults = struct( ...
        "well",    num2cell(wells), ...
        "Neurons", num2cell(Neurons, 2), ...
        "Glia",    num2cell(Glia, 2) ...
        );
    % For full per-cell stats across wells, you can also store
    % Neuron_stats_all / Glia_stats_all in a similar fashion.

    distances = struct( ...
        "NtoN",  AllNtoN,  ...
        "GtoG",  AllGtoG,  ...
        "NtoC",  AllNtoC,  ...
        "GtoC",  AllGtoC,  ...
        "NtoNC", AllNtoNC, ...
        "GtoGC", AllGtoGC, ...
        "NtoNL", AllNtoNL, ...
        "GtoGL", AllGtoGL, ...
        "NtoGL", AllNtoGL, ...
        "GtoNL", AllGtoNL  ...
        );

    cellResultData(i, :) = {dataset, datasetResults, cellResultWells, distances};

    %% Save dataset-level results to disk

    datasetSaveDir = fullfile(savepath, dataset);
    if ~isfolder(datasetSaveDir)
        mkdir(datasetSaveDir);
    end

    parsave(fullfile(datasetSaveDir, [char(dataset) '.mat']), datasetResults);
    parsave(fullfile(datasetSaveDir, 'distances.mat'), distances);

    disp("Completed dataset: " + dataset);

end % datasets parfor

toc;
delete(gcp('nocreate'));
disp("Finished all datasets.");
disp("Script complete!");
