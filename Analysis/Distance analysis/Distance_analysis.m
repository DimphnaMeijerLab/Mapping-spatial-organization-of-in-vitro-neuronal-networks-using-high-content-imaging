function [NtoN, GtoG, NtoC, GtoC, NtoNC, GtoGC, NtoNL, GtoGL, NtoGL, GtoNL] = Distance_analysis(neuronCentroids, gliaCentroids)
% Distance_analysis: Computes distances between neurons and glial cells.
% 
% Inputs:
%   neuronCentroids - Nx2 matrix of neuron centroid coordinates (x, y)
%   gliaCentroids   - Mx2 matrix of glia centroid coordinates (x, y)
% 
% Outputs:
%   NtoN  - Distances between all neuron pairs
%   GtoG  - Distances between all glia pairs
%   NtoC  - Distances from neurons to the image center
%   GtoC  - Distances from glia to the image center
%   NtoNC - Neuron-to-neuron distances within half image radius
%   GtoGC - Glia-to-glia distances within half image radius
%   NtoNL - Neuron-to-neuron distances, excluding those near the edge
%   GtoGL - Glia-to-glia distances, excluding those near the edge
%   NtoGL - Neuron-to-glia distances
%   GtoNL - Glia-to-neuron distances

%% Constants
image_radius = 7750; % Half of the image width (assuming square image)
pixel_to_um = 878 / 1104; % Corrected conversion factor from pixels to micrometers
edge_buffer = 700 * pixel_to_um; % Edge exclusion distance in micrometers

%% Identify central and non-central cells
half_radius = image_radius / 2; % Define half-radius threshold
central_neurons = neuronCentroids(vecnorm(neuronCentroids, 2, 2) < half_radius, :);
central_glia = gliaCentroids(vecnorm(gliaCentroids, 2, 2) < half_radius, :);
non_central_neurons = neuronCentroids(vecnorm(neuronCentroids, 2, 2) >= half_radius, :);
non_central_glia = gliaCentroids(vecnorm(gliaCentroids, 2, 2) >= half_radius, :);

%% Identify local and non-local cells
local_radius = image_radius - edge_buffer; % Define local cell radius
local_neurons = neuronCentroids(vecnorm(neuronCentroids, 2, 2) < local_radius, :);
local_glia = gliaCentroids(vecnorm(gliaCentroids, 2, 2) < local_radius, :);
non_local_neurons = neuronCentroids(vecnorm(neuronCentroids, 2, 2) >= local_radius, :);
non_local_glia = gliaCentroids(vecnorm(gliaCentroids, 2, 2) >= local_radius, :);

%% Compute distances and convert to micrometers
NtoN = pdist(neuronCentroids) * pixel_to_um;
GtoG = pdist(gliaCentroids) * pixel_to_um;
NtoC = vecnorm(neuronCentroids, 2, 2)' * pixel_to_um;
GtoC = vecnorm(gliaCentroids, 2, 2)' * pixel_to_um;

NtoNC = [reshape(pdist2(central_neurons, non_central_neurons).', 1, []), pdist(central_neurons)] * pixel_to_um;
GtoGC = [reshape(pdist2(central_glia, non_central_glia).', 1, []), pdist(central_glia)] * pixel_to_um;

NtoNL = [reshape(pdist2(local_neurons, non_local_neurons).', 1, []), pdist(local_neurons)] * pixel_to_um;
GtoGL = [reshape(pdist2(local_glia, non_local_glia).', 1, []), pdist(local_glia)] * pixel_to_um;

NtoGL = pdist2(neuronCentroids, gliaCentroids) * pixel_to_um;
GtoNL = pdist2(gliaCentroids, neuronCentroids) * pixel_to_um;

%% Remove self-distances (zero values)
NtoN(NtoN == 0) = [];
GtoG(GtoG == 0) = [];
NtoC(NtoC == 0) = [];
GtoC(GtoC == 0) = [];
NtoNC(NtoNC == 0) = [];
GtoGC(GtoGC == 0) = [];
NtoNL(NtoNL == 0) = [];
GtoGL(GtoGL == 0) = [];

%% Remove excessively long distances
NtoNC(NtoNC >= half_radius * pixel_to_um) = [];
GtoGC(GtoGC >= half_radius * pixel_to_um) = [];
NtoNL(NtoNL >= 500) = [];
GtoGL(GtoGL >= 100) = [];
NtoGL(NtoGL >= 100) = 0;
GtoNL(GtoNL >= 100) = 0;

end