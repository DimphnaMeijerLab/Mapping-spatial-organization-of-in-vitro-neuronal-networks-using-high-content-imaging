%% Thank you for using Mia
% This is Mia, the Modular Image Analysis tool that can analyse 
% large datasets of whole well images simultaneously using parfor!
% To find out how to set up your file structure, please check the 
% Tutorial on GitHub
%
% Look for the designated
%% Comment sections 
% to find out what you need to do to add your own modules and to help you
% get set up

%init
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

%% Data selection

%Input path
disp('Select the folder that contains your Data')
datapath = uigetdir('E:\', ...
    'Select the folder that contains your Data'); 

%List Datasets in Data folder
tempdatasets = GetSubDirsFirstLevelOnly(datapath);

%Select desired Folders/Datasets in Data folder
include = listdlg( ...
    'ListString', ...
        tempdatasets, ...
    'PromptString', ...
        'Please select the folder(s) that you want to analyse', ...
    'ListSize', ...
        [300,300] ...
    );

%Select Datasets
datasets = [];
for i = 1:length(tempdatasets)
    dataset = string(tempdatasets(i));
    if ismember(dataset, tempdatasets(include))
        datasets = [datasets; dataset]; %#ok<AGROW> 
    end
end

%Output path
disp("Select the folder in which you'd like to save your Results")
savepath = uigetdir('E:\', ...
    "Select the folder in which you'd like to save your Results");


%% Enter Worker count
workercount = 16;

%----Parfor loop over Datasets----
tic
parpool('local',workercount);
cellResultData = cell(length(datasets),4);

parfor i = 1:length(datasets) %parfor

    dataset = string(datasets(i));
    disp("Doing Dataset: " + dataset);

    %list Wells, make string from cell array and transpose for data struct
    wells = string(GetSubDirsFirstLevelOnly(datapath + "\" + dataset))';

    % Initialize module variables
    % We initialise with the correct datatype for a slight speedup

    LBW_Dapi = false(2); 
    LBW_AnkG = false(2); 
    LBW_Neurons = false(2);
    
        % For nuclei analysis
    Neurons = zeros(length(wells),1);
    Neuron_stats = struct();
    Neuron_stats_all = cell(1,length(wells));
    Glia = zeros(length(wells),1);
    Glia_stats = struct();
    Glia_stats_all = cell(1,length(wells));
    Thresholds = zeros(length(wells),1);

        % For distance analysis
    NtoN = zeros(1);
    GtoG = zeros(1);
    NtoC = zeros(1);
    GtoC = zeros(1);
    NtoNC = zeros(1);
    GtoGC = zeros(1);
    NtoNL = zeros(1);
    GtoGL = zeros(1);
    NtoGL = zeros(1);
    GtoNL = zeros(1);
    AllNtoN = zeros(1);
    AllGtoG = zeros(1);
    AllNtoC = zeros(1);
    AllGtoC = zeros(1);
    AllNtoNC = zeros(1);
    AllGtoGC = zeros(1);
    AllNtoNL = zeros(1);
    AllGtoGL = zeros(1);
    AllNtoGL = zeros(1);
    AllGtoNL = zeros(1);

    cellResultWells = cell(length(wells),2);
    %----Loop over Wells----
    for j = 1:length(wells)

        well = wells(j);
        disp('Doing well: ' + well + ' of ' + dataset);

        % Load Required Image Files
        
        %Dapi analysis
        LBW_Dapi = logical(imread(strcat( ...
            datapath,'\',dataset,'\',well,'\',well,'_th_dapi.tif')));
        LBW_AnkG = logical(imread(strcat( ...
            datapath,'\',dataset,'\',well,'\',well,'_th_AnkG.tif')));
        Dapi = imread(strcat( ...
            datapath,'\',dataset,'\',well,'\',well,'_dapi.tif'));

        % Add Modules
        % Here you can add your modules, make sure your script outputs
        % the variables you care about as a list/Mx1 dimensional matrix
        % using [ ]

        % Example: 
        %           [   
        %               Glia_stats, Neuron_stats, ...
        %               LBW_Glia, LBW_Neurons, ...
        %               Thresholds(j,1) ...
        %           ]   = Nuclei_analysis(LBW_Dapi,LBW_AnkG,Dapi);
        %
        % Don't forget to put your output in a variable that can store the
        % data for the entire multiplicate(dataset) properly

        % Example: 
        %             Neurons(j)=length(Neuron_stats(:,1));
        %             Glia(j)=length(Glia_stats(:,1));
        %


        % Module 1 - Nuclei Analysis
        disp("Analysing Nuclei...");
        [   
            Glia_stats, Neuron_stats, ...
            LBW_Glia, LBW_Neurons, ...
            Thresholds(j,1) ...
        ]   = Nuclei_analysis(LBW_Dapi, LBW_AnkG,Dapi);

        %Count Neurons and Glia
        Neurons(j) = length(Neuron_stats.Area);
        Glia(j) = length(Glia_stats.Area);
        
%         disp(Neurons(j));
%         disp(Glia(j));

        disp("Nuclei Analysis complete");
        disp("Calculating distances...");

        % Module 2 - Distance Analysis
        [
            NtoN, GtoG, NtoC, ...
            GtoC, NtoNC, GtoGC, ...
            NtoNL, GtoGL, NtoGL, GtoNL ...
        ] = Distance_analysis(Neuron_stats.Centroid, Glia_stats.Centroid);

        AllNtoN = [AllNtoN, NtoN];
        AllGtoG = [AllGtoG, GtoG];
        AllNtoC = [AllNtoC, NtoC];
        AllGtoC = [AllGtoC, GtoC];
        AllNtoNC = [AllNtoNC, NtoNC];
        AllGtoGC = [AllGtoGC, GtoGC];
        AllNtoNL = [AllNtoNL, NtoNL];
        AllGtoGL = [AllGtoGL, GtoGL];
%         AllNtoGL = [AllNtoGL, NtoGL];
%         AllGtoNL = [AllGtoNL, GtoNL];

        disp("Distance Analysis complete");
    
        disp("Analysis complete of well " + well + " of " + dataset);


        % Save well results
        % Here you save your well results in a struct and name them.

        % Example:
        %       wellResults= ...
        %       struct  ( ...
        %               "Neurons", LBW_Neurons, ...
        %               "Glia", LBW_Glia, ...
        %               "Neuronstats", Neuron_stats, ...
        %               "Gliastats", Glia_stats ...
        %               );

        wellResults = ...
        struct  ( ...
                "Neurons", LBW_Neurons, ...
                "Glia", LBW_Glia, ...
                "Neuronstats", Neuron_stats, ...
                "Gliastats", Glia_stats, ...
                'NtoN', NtoN, ...
                'GtoG', GtoG, ...
                'NtoC', NtoC, ...
                'GtoC', GtoC, ...
                'NtoNC', NtoNC, ...
                'GtoGC', GtoGC, ...
                'NtoNL', NtoNL, ...
                'GtoGL', GtoGL, ...
                'NtoGL', NtoGL, ...
                'GtoNL', GtoNL ...
                );
        
%         if j == 1
%             Neuron_stats_all = Neuron_stats;
%             Glia_stats_all = Glia_stats;
%         else
%             Neuron_stats_all = [Neuron_stats_all Neuron_stats];
%             Glia_stats_all = [Glia_stats_all Glia_stats];
%         end

        cellResultWells(j,:) = {well, wellResults};

        % Saving well results
        % First check if the desired directory is there, and create
        % it if it is not
        if not(isfolder(strcat(savepath,'\',dataset,'\',well,'\')))
            mkdir(strcat(savepath,'\',dataset,'\',well,'\'));
        end
        % Save the well results
        parsave(strcat(savepath,'\',dataset,'\',well,'\', ...
            well,'.mat'),wellResults);

    end
    
    % Save Dataset Results
    % Here you save your Dataset results in a struct and name them.
    % Each element in the struct needs to be a cell and the first element 
    % should be the list of wells. Don't forget that the final element
    % doesnt need a comma at the end.

    % Example:
    %       datasetResults = ...
    %       struct  ( ...
    %               "well", num2cell(wells), ...
    %               "Neurons", num2cell(Neurons,2), ...
    %               "Glia", num2cell(Glia,2) ...
    %               );

    datasetResults = ...
    struct  ( ...
            "well", num2cell(wells), ...
            "Neurons", num2cell(Neurons,2), ...
            "Glia", num2cell(Glia,2) ...
            );
%             "Neuron_stats", num2cell(Neuron_stats_all)', ...
%             "Glia_stats", num2cell(Glia_stats_all)' ...

    distances = ...
    struct  ( ...
            "NtoN", AllNtoN, ...
            "GtoG", AllGtoG, ...
            "NtoC", AllNtoC, ...
            "GtoC", AllGtoC, ...
            "NtoNC", AllNtoNC, ...
            "GtoGC", AllGtoGC, ...
            "NtoNL", AllNtoNL, ...
            "GtoGL", AllGtoGL, ...
            "NtoGL", AllNtoGL, ...
            "GtoNL", AllGtoNL ...
            );

    cellResultData(i,:) = {dataset, datasetResults, cellResultWells, distances};

    % Saving well results
    % First check if the desired directory is there, and create
    % it if it is not
    if not(isfolder(strcat(savepath,'\',dataset,'\')))
        mkdir(strcat(savepath,'\',dataset,'\'));
    end
    % Save the well results
    parsave(strcat(savepath,'\',dataset,'\', ...
        dataset,'.mat'),datasetResults);
    parsave(strcat(savepath,'\',dataset,'\', ...
        'distances.mat'),distances);
    
    disp("Done with " + dataset)
end

toc
delete(gcp('nocreate'));
disp("Finished all datasets")

% %% ----save parfor----
% 
% tic
% disp("Now saving all dataset results...");
% 
% datasets = [cellResultData{:,1}];
% for i=1:length(datasets)
%     dataset = cellResultData{i,1};
%     datasetResults = cellResultData{i,2};
%     if not(isfolder(strcat(savepath,'\',dataset,'\')))
%         mkdir(strcat(savepath,'\',dataset,'\'));
%     end
%     save(strcat(savepath,'\',dataset,'\', ...
%         dataset,'.mat'),'datasetResults','-v7.3');
% 
%     cellResultWells = cellResultData{i,3};
%     wells=[cellResultWells{:,1}];
%     for j=1:length(wells)
%         well = cellResultWells{j,1};
%         wellResults = cellResultWells{j,2};
%         if not(isfolder(strcat(savepath,'\',dataset,'\',well,'\')))
%             mkdir(strcat(savepath,'\',dataset,'\',well,'\'));
%         end
%         save(strcat(savepath,'\',dataset,'\',well,'\', ...
%             well,'.mat'),'wellResults','-v7.3');
%     end
% end
% disp("All dataset results saved.");
disp("Script Complete!");
% toc

