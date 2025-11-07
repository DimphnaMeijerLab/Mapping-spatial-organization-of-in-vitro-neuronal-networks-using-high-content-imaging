function [NtoN, GtoG, NtoC, GtoC, NtoNC, GtoGC, NtoNL, GtoGL, NtoGL, GtoNL] = ...
    Distance_analysis(neuronCentroids, gliaCentroids, pixelSize, wellRadius)
%DISTANCE_ANALYSIS Compute distance metrics between neuronal and glial nuclei.
%
%   [NtoN, GtoG, NtoC, GtoC, NtoNC, GtoGC, NtoNL, GtoGL, NtoGL, GtoNL] = ...
%       DISTANCE_ANALYSIS(neuronCentroids, gliaCentroids, pixelSize, wellRadius)
%
% Inputs:
%   neuronCentroids - Nx2 matrix of neuron centroid coordinates (x, y) in pixels
%   gliaCentroids   - Mx2 matrix of glia centroid coordinates (x, y) in pixels
%   pixelSize       - Scalar, linear pixel size in µm/pixel
%   wellRadius      - Scalar, physical radius of the well (or image) in µm
%
% Outputs (all distances in µm):
%   NtoN  - Distances between all neuron pairs
%   GtoG  - Distances between all glia pairs
%   NtoC  - Distances from neurons to the well/image center
%   GtoC  - Distances from glia to the well/image center
%   NtoNC - Neuron-to-neuron distances within half well radius
%   GtoGC - Glia-to-glia distances within half well radius
%   NtoNL - Neuron-to-neuron distances, excluding nuclei near the edge
%   GtoGL - Glia-to-glia distances, excluding nuclei near the edge
%   NtoGL - Neuron-to-glia distances
%   GtoNL - Glia-to-neuron distances
%
% Assumptions:
%   - Centroid coordinates are expressed relative to the well/image center,
%     i.e. (0,0) is center. If they are in image coordinates (origin top-left),
%     recenter before calling this function.
%
% All thresholds and radii are applied in physical units (µm).

%% Convert input to double and physical units (µm)

neuronCentroids = double(neuronCentroids);
gliaCentroids   = double(gliaCentroids);

% Convert from pixels to µm
neuronUm = neuronCentroids * pixelSize;
gliaUm   = gliaCentroids   * pixelSize;

% Radial distances from center (µm)
neuronRadiusUm = vecnorm(neuronUm, 2, 2);
gliaRadiusUm   = vecnorm(gliaUm,   2, 2);

%% Define radii and thresholds (all in µm)

halfRadiusUm = wellRadius / 2;          % "central" region radius
edgeBufferUm = 700 * pixelSize;         % ~700 pixels from well edge, in µm
localRadiusUm = max(wellRadius - edgeBufferUm, 0);

% Additional cut-offs for long distances (µm)
maxNNLocalUm = 500;   % for NtoNL
maxGGLocalUm = 100;   % for GtoGL
maxNGUm      = 100;   % for NtoGL
maxGNUm      = 100;   % for GtoNL

%% Identify central and non-central cells (by half well radius)

central_neurons     = neuronUm(neuronRadiusUm <  halfRadiusUm, :);
non_central_neurons = neuronUm(neuronRadiusUm >= halfRadiusUm, :);

central_glia        = gliaUm(gliaRadiusUm <  halfRadiusUm, :);
non_central_glia    = gliaUm(gliaRadiusUm >= halfRadiusUm, :);

%% Identify local and non-local cells (exclude edge region)

local_neurons     = neuronUm(neuronRadiusUm <  localRadiusUm, :);
non_local_neurons = neuronUm(neuronRadiusUm >= localRadiusUm, :);

local_glia        = gliaUm(gliaRadiusUm <  localRadiusUm, :);
non_local_glia    = gliaUm(gliaRadiusUm >= localRadiusUm, :);

%% Compute distances (all in µm)

% Neuron-neuron and glia-glia
NtoN = pdist(neuronUm);   % row vector
GtoG = pdist(gliaUm);     % row vector

% Distances to center
NtoC = neuronRadiusUm.';  % row vector
GtoC = gliaRadiusUm.';    % row vector

% Central vs non-central + within-central (neuron)
distCentralNonCentralN = pdist2(central_neurons, non_central_neurons);
NtoNC = [ ...
    reshape(distCentralNonCentralN.', 1, []), ...
    pdist(central_neurons) ...
    ];

% Central vs non-central + within-central (glia)
distCentralNonCentralG = pdist2(central_glia, non_central_glia);
GtoGC = [ ...
    reshape(distCentralNonCentralG.', 1, []), ...
    pdist(central_glia) ...
    ];

% Local vs non-local + within-local (neuron)
distLocalNonLocalN = pdist2(local_neurons, non_local_neurons);
NtoNL = [ ...
    reshape(distLocalNonLocalN.', 1, []), ...
    pdist(local_neurons) ...
    ];

% Local vs non-local + within-local (glia)
distLocalNonLocalG = pdist2(local_glia, non_local_glia);
GtoGL = [ ...
    reshape(distLocalNonLocalG.', 1, []), ...
    pdist(local_glia) ...
    ];

% Neuron-glia and glia-neuron
NtoGL = pdist2(neuronUm, gliaUm);   % matrix (neurons x glia)
GtoNL = pdist2(gliaUm,   neuronUm); % matrix (glia x neurons)

% Vectorise NtoGL and GtoNL for further analysis
NtoGL = NtoGL(:).';  % row vector
GtoNL = GtoNL(:).';  % row vector

%% Remove exact zeros (self distances / degenerate cases)

epsilon = 1e-9;

NtoN(abs(NtoN) < epsilon) = [];
GtoG(abs(GtoG) < epsilon) = [];
NtoC(abs(NtoC) < epsilon) = [];
GtoC(abs(GtoC) < epsilon) = [];
NtoNC(abs(NtoNC) < epsilon) = [];
GtoGC(abs(GtoGC) < epsilon) = [];
NtoNL(abs(NtoNL) < epsilon) = [];
GtoGL(abs(GtoGL) < epsilon) = [];
NtoGL(abs(NtoGL) < epsilon) = [];
GtoNL(abs(GtoNL) < epsilon) = [];

%% Remove excessively long distances (cut-offs in µm)

% Central vs non-central: cut at half well radius
NtoNC(NtoNC >= halfRadiusUm) = [];
GtoGC(GtoGC >= halfRadiusUm) = [];

% Local vs non-local (neuron/glia)
NtoNL(NtoNL >= maxNNLocalUm) = [];
GtoGL(GtoGL >= maxGGLocalUm) = [];

% Neuron-glia cross distances
NtoGL(NtoGL >= maxNGUm) = [];
GtoNL(GtoNL >= maxGNUm) = [];

end
