function cluster_stats = cluster_analysis(neuron_centroids, critical_distance_um)
% CLUSTER_ANALYSIS
%   Analyze clusters of neurons based on their centroids (in pixels).
%
%   cluster_stats = cluster_analysis(neuron_centroids, critical_distance_um)
%
%   Input:
%       neuron_centroids      - Nx2 matrix of neuron centroids ([x, y] in pixels)
%       critical_distance_um  - distance threshold (µm) to consider neurons neighbors
%
%   Output:
%       cluster_stats - struct with:
%           .num_clusters             - number of connected components (incl. singletons)
%           .cluster_sizes            - vector of cluster sizes (one entry per component)
%           .num_singletons           - number of neurons in size-1 components
%           .num_clusters_gt1         - number of components with size > 1
%           .mean_cluster_size_gt1    - mean size of clusters > 1
%           .median_cluster_size_gt1  - median size of clusters > 1
%           .fraction_in_clusters_gt1 - fraction of neurons in clusters > 1
%           .cluster_ids              - component index for each neuron
%           .critical_distance_um     - threshold used (µm)

    % --- Conversion factor µm -> pixels (adapt to your imaging setup if needed) ---
    um_to_pixel = 1104 / 878;   % your original scale
    critical_distance_px = critical_distance_um * um_to_pixel;

    % --- Basic checks ---
    if isempty(neuron_centroids) || size(neuron_centroids, 2) ~= 2
        error('neuron_centroids must be an Nx2 matrix of [x, y] coordinates.');
    end

    num_neurons = size(neuron_centroids, 1);

    % --- Pairwise distances in pixels ---
    dist_matrix = pdist2(neuron_centroids, neuron_centroids);

    % --- Adjacency: 1 if distance < threshold (exclude self) ---
    adjacency_matrix = sparse(dist_matrix < critical_distance_px & dist_matrix > 0);

    % --- Graph & connected components ---
    G = graph(adjacency_matrix);

    % Correct usage of conncomp:
    %  - first output: component index for each node
    %  - second output: size of each component
    [cluster_ids_row, component_sizes] = conncomp(G);
    cluster_ids   = cluster_ids_row(:);    % Nx1
    cluster_sizes = component_sizes(:);    % Cx1, C = number of components

    % Number of clusters = number of components (incl. singletons)
    num_clusters = numel(cluster_sizes);

    % --- Stats ---
    num_singletons = sum(cluster_sizes == 1);
    num_in_clusters_gt1 = num_neurons - num_singletons;
    fraction_in_clusters_gt1 = num_in_clusters_gt1 / num_neurons;

    clusters_gt1 = cluster_sizes(cluster_sizes > 1);
    if isempty(clusters_gt1)
        num_clusters_gt1        = 0;
        mean_cluster_size_gt1   = NaN;
        median_cluster_size_gt1 = NaN;
    else
        num_clusters_gt1        = numel(clusters_gt1);
        mean_cluster_size_gt1   = mean(clusters_gt1);
        median_cluster_size_gt1 = median(clusters_gt1);
    end

    % --- Output struct ---
    cluster_stats = struct( ...
        'num_clusters',             num_clusters, ...
        'cluster_sizes',            cluster_sizes, ...
        'num_singletons',           num_singletons, ...
        'num_clusters_gt1',         num_clusters_gt1, ...
        'mean_cluster_size_gt1',    mean_cluster_size_gt1, ...
        'median_cluster_size_gt1',  median_cluster_size_gt1, ...
        'fraction_in_clusters_gt1', fraction_in_clusters_gt1, ...
        'cluster_ids',              cluster_ids, ...
        'critical_distance_um',     critical_distance_um ...
    );
end

