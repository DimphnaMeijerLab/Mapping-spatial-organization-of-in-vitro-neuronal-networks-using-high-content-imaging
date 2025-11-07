% =========================================================================
% Main_Cluster_analysis_v1.m
% -------------------------------------------------------------------------
% - Recursively scans a root folder for all .mat files
% - For each file:
%       * checks there is a valid well ID in the path
%       * only keeps wells in the set {B02, C02, D02, B08, C08, D08}
%       * loads x.Neuronstats.Centroid (Nx2, in pixels)
%       * runs cluster_analysis() for that file
% - Extracts DIV from the path
% - Treats singletons as isolated neurons, not real clusters
% - Stores per file:
%       * DIV, WellID, FileName, NumNeurons
%       * TotalClusters (incl. singletons)
%       * NumClusters_gt1 (clusters with size > 1)
%       * NumSingletons
%       * MeanClusterSize_gt1, MedianClusterSize_gt1
%       * FractionInClusters_gt1
%       * ClusterSizes (full vector)
% =========================================================================

clear; clc; close all;

%% === 1. Settings ===

% Root folder containing all DIVs, pups, wells, etc.
root_dir = 'E:\Angelica\Results\CerebellarNeurons\AC011_DEFINITIVA\500um';  % <-- adjust if needed

% Critical distance (µm) for clustering
critical_distance_um = 110;   % <-- adjust based on your NN-distance peak

% Only analyze these wells:
allowed_wells = {'B04','C03','B03','C04'};

%% === 2. Find all .mat files recursively ===

search_pattern = fullfile(root_dir, '**', '*.mat');
files = dir(search_pattern);

if isempty(files)
    error('No .mat files found under %s', root_dir);
end

fprintf('Found %d .mat files.\n\n', numel(files));

%% === 3. Loop over files and run cluster analysis PER FILE ===

file_names        = {};
full_paths        = {};
div_labels        = {};
well_labels       = {};
num_neurons         = [];
total_clusters    = [];  % including singletons
num_clusters_gt1  = [];  % only clusters with size > 1
num_singletons    = [];
mean_clust_size_gt1   = [];
median_clust_size_gt1 = [];
fraction_gt1      = [];
cluster_sizes_all = {};  % full cluster size vector per file

for iFile = 1:numel(files)
    file_path = fullfile(files(iFile).folder, files(iFile).name);

    % --- Detect DIV from full path (case-insensitive: Div7, DIV7, div_7, etc.) ---
    tokens_div = regexp(file_path, '(?i)div[_-]?(\d+)', 'tokens', 'once');
    if isempty(tokens_div)
        div_label = 'Unknown_DIV';
    else
        div_label = ['Div' tokens_div{1}];   % e.g. 'Div7'
    end

    % --- Detect WELL ID from path (A01..H12, etc.) ---
    % Look for "\A01\" or "/A01/" style segments
    tokens_well = regexp(file_path, '[\\/](?<well>[A-H][0-9]{2})[\\/]', 'names', 'once');
    if isempty(tokens_well)
        fprintf('Skipping (no well ID like A01/B02 found in path): %s\n', file_path);
        continue;
    end
    well_label = upper(tokens_well.well);  % e.g. 'B02'

    % --- Keep only selected wells ---
    if ~ismember(well_label, allowed_wells)
        fprintf('Skipping well %s (not in allowed set).\n', well_label);
        continue;
    end

    % --- Load data ---
    data = load(file_path);

    %Expect x.neuronstats.Centroid as Nx2 (in pixels)
    if ~isfield(data, 'var') || ~isfield(data.var, 'NeuronCentroids') 
       %~isfield(data.x.Neuronstats, 'Centroid')
        %fprintf('Skipping (no x.Neuronstats.Centroid): %s\n', file_path);
        continue;
    end

    neuron_centroids = data.var.NeuronCentroids;

    if size(neuron_centroids, 2) ~= 2
        fprintf('Skipping (Centroid not Nx2): %s\n', file_path);
        continue;
    end

    num_neuron_in_file = size(neuron_centroids, 1);

    % --- Run cluster analysis for THIS file ---
    stats = cluster_analysis(neuron_centroids, critical_distance_um);

    % Pull out relevant fields
    this_total_clusters      = stats.num_clusters;              % incl. singletons
    this_num_clusters_gt1    = stats.num_clusters_gt1;          % only clusters > 1
    this_num_singletons      = stats.num_singletons;
    this_mean_cs_gt1         = stats.mean_cluster_size_gt1;
    this_median_cs_gt1       = stats.median_cluster_size_gt1;
    this_fraction_gt1        = stats.fraction_in_clusters_gt1;  % neurons in clusters > 1
    this_cluster_sizes       = stats.cluster_sizes;             % full vector

    % --- Store in arrays / cell arrays ---
    file_names{end+1,1}        = files(iFile).name;
    full_paths{end+1,1}        = file_path;
    div_labels{end+1,1}        = div_label;
    well_labels{end+1,1}       = well_label;

    num_neurons(end+1,1)       = num_neuron_in_file;
    total_clusters(end+1,1)    = this_total_clusters;
    num_clusters_gt1(end+1,1)  = this_num_clusters_gt1;
    num_singletons(end+1,1)    = this_num_singletons;
    mean_clust_size_gt1(end+1,1)   = this_mean_cs_gt1;
    median_clust_size_gt1(end+1,1) = this_median_cs_gt1;
    fraction_gt1(end+1,1)      = this_fraction_gt1;
    cluster_sizes_all{end+1,1} = this_cluster_sizes;

    % --- Inline print (readable summary) ---
    fprintf(['% -8s | Well: %-4s | %-30s | Neurons: %4d | Total cl: %4d | ' ...
             'Cl>1: %4d | mean size>1: %.2f | median size>1: %.2f | ' ...
             'frac in cl>1: %.2f | singletons: %4d\n'], ...
        div_label, well_label, files(iFile).name, num_neuron_in_file, ...
        this_total_clusters, this_num_clusters_gt1, ...
        this_mean_cs_gt1, this_median_cs_gt1, ...
        this_fraction_gt1, this_num_singletons);
end

%% === 4. Combine results into table and sort ===

if isempty(file_names)
    error('No valid centroid files processed. Check paths, wells and x.Neuronstats.Centroid.');
end

results = table( ...
    div_labels, well_labels, file_names, full_paths, ...
    num_neurons, total_clusters, num_clusters_gt1, num_singletons, ...
    mean_clust_size_gt1, median_clust_size_gt1, fraction_gt1, ...
    cluster_sizes_all, ...
    'VariableNames', {'DIV', 'WellID', 'FileName', 'FullPath', ...
                      'NumNeurons', 'TotalClusters', 'NumClusters_gt1', 'NumSingletons', ...
                      'MeanClusterSize_gt1', 'MedianClusterSize_gt1', 'FractionInClusters_gt1', ...
                      'ClusterSizes'} );

% Sort by DIV, then Well, then FileName
results = sortrows(results, {'DIV', 'WellID', 'FileName'});

fprintf('\n=== Full results table (sorted by DIV, Well, FileName) ===\n');
disp(results);

%% === 5. Grouped summary by DIV and Well ===

unique_divs = unique(results.DIV, 'stable');
fprintf('\n=== Grouped summary by DIV and Well ===\n');

for iDiv = 1:numel(unique_divs)
    this_div = unique_divs{iDiv};
    idx_div = strcmp(results.DIV, this_div);

    fprintf('\n----- %s -----\n', this_div);

    sub_div = results(idx_div, :);

    unique_wells = unique(sub_div.WellID, 'stable');
    for iWell = 1:numel(unique_wells)
        this_well = unique_wells{iWell};
        idx_well = strcmp(sub_div.WellID, this_well);
        sub = sub_div(idx_well, :);

        fprintf('  Well %s:\n', this_well);
        for iRow = 1:height(sub)
            fprintf(['    File: %-30s | Neurons: %4d | Total cl: %4d | ' ...
                     'Cl>1: %4d | mean size>1: %.2f | median size>1: %.2f | ' ...
                     'frac in cl>1: %.2f | singletons: %4d\n'], ...
                sub.FileName{iRow}, sub.NumNeurons(iRow), sub.TotalClusters(iRow), ...
                sub.NumClusters_gt1(iRow), sub.MeanClusterSize_gt1(iRow), ...
                sub.MedianClusterSize_gt1(iRow), sub.FractionInClusters_gt1(iRow), ...
                sub.NumSingletons(iRow));
        end
    end
end

%% === 6. Save results ===

out_csv = fullfile(root_dir, 'cluster_counts_and_sizes_by_DIV_and_Well_selectedWells.csv');
writetable(results, out_csv);
fprintf('\nSaved detailed results to: %s\n', out_csv);

out_mat = fullfile(root_dir, 'cluster_counts_and_sizes_by_DIV_and_Well_selectedWells.mat');
save(out_mat, 'results');
fprintf('Saved MATLAB results to: %s\n', out_mat);

