function [NtoN, NtoC, NtoNC, NtoNL] = Distance_analysis_centroids_v3(neuronCentroids)
    % Calculate distances between neurons using their centroids
    
    % Image parameters
    imageDiameter = 18976; % Example diameter in pixels, modify if different
    imageRadius = imageDiameter / 2; % Radius of the image in pixels
    pixelToMicrometerRatio = 1.257; % Conversion ratio from pixels to micrometers
    micrometersToPixels = 1 / pixelToMicrometerRatio; % Conversion from micrometers to pixels
    
    % Define the actual well center coordinates (e.g., [x, y])
    wellCenter = [9488, 9488]; % Modify this to the actual center of the well in pixels
    
    % Distance threshold in micrometers (converted to pixels)
    distanceThresholdMicrometers = 100; % Distance from the edge in micrometers
    distanceThresholdPixels = distanceThresholdMicrometers * micrometersToPixels; % Convert to pixels
    
    % Calculate distances from each neuron to the well center
    neuronDistancesToWellCenter = sqrt(sum((neuronCentroids - wellCenter).^2, 2));
    
    % Find central neurons within half the well radius
    central_neurons_locations = neuronDistancesToWellCenter < (imageRadius / 2);
    central_neurons = neuronCentroids(central_neurons_locations, :);
    
    % Find local neurons not closer than a specific distance to the edge
    local_neurons_locations = neuronDistancesToWellCenter < (imageRadius - distanceThresholdPixels);
    local_neurons = neuronCentroids(local_neurons_locations, :);
    
    % Find non-central neurons
    non_central_neurons = neuronCentroids(~central_neurons_locations, :);
    
    % Find non-local neurons
    non_local_neurons = neuronCentroids(~local_neurons_locations, :);
    
    % Calculating the distances and converting from pixels to micrometers
    NtoN = pdist(neuronCentroids) * micrometersToPixels;
    NtoC = neuronDistancesToWellCenter' * micrometersToPixels;
    
    if size(central_neurons, 1) > 1
        NtoNC = [reshape( ...
            pdist2(central_neurons, non_central_neurons).',1,[]), ...
            pdist(central_neurons)] * micrometersToPixels;
    else
        NtoNC = []; % Not enough central neurons to calculate distances
    end
    
    if size(local_neurons, 1) > 1
        NtoNL = [reshape( ...
            pdist2(local_neurons, non_local_neurons).',1,[]), ...
            pdist(local_neurons)] * micrometersToPixels;
    else
        NtoNL = []; % Not enough local neurons to calculate distances
    end

    % Removing distance to itself (==0)
    NtoN(NtoN == 0) = [];
    NtoC(NtoC == 0) = [];
    NtoNC(NtoNC == 0) = [];
    NtoNL(NtoNL == 0) = [];

    % Do not need to save distances that are too long
    NtoNC(NtoNC >= imageRadius) = [];
    NtoNL(NtoNL >= distanceThresholdMicrometers) = [];
end