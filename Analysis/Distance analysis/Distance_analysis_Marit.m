function [NtoN, GtoG, NtoC, GtoC, NtoNC, GtoGC, NtoNL, GtoGL, NtoGL, GtoNL] = Distance_analysis_Marit(neuronCentroids, gliaCentroids)
% This function calculates distances between neurons and glia cells
    % Using the centroids (where the mid point of the image is (0,0)) of 
    % each neuron/glia cell, this function will calculate the distance 
    % between: 
        % NtoN: neurons to neurons
        % GtoG: glia to glia
        % NtoC: neurons to the center of the image
        % GtoC: glia to the center of the image
        % NtoNC: neurons to neurons, but only within half the radius of
            % the image
        % GtoGC: glia to glia, but only within half the radius of the 
            % image
    
    % Find the cells with their centroid within half the radius
    central_neurons_locations = sqrt(neuronCentroids(:,1).^2 + ...
        neuronCentroids(:,2).^2) < (12006/4);
    central_glia_locations = sqrt(gliaCentroids(:,1).^2 + ...
        gliaCentroids(:,2).^2) < (12006/4);
    central_neurons = neuronCentroids(central_neurons_locations,:);
    central_glia = gliaCentroids(central_glia_locations,:);
    non_central_neurons = neuronCentroids(not(central_neurons_locations),:);
    non_central_glia = gliaCentroids(not(central_glia_locations),:);

    % Find the cells with their centroid not closer than 100 um to the edge
    local_neurons_locations = sqrt(neuronCentroids(:,1).^2 + ...
        neuronCentroids(:,2).^2) < (12006/2-(700*12006/4620));
    local_glia_locations = sqrt(gliaCentroids(:,1).^2 + ...
        gliaCentroids(:,2).^2) < (12006/2-(700*12006/4620));
    local_neurons = neuronCentroids(local_neurons_locations,:);
    local_glia = gliaCentroids(local_glia_locations,:);
    non_local_neurons = neuronCentroids(not(local_neurons_locations),:);
    non_local_glia = gliaCentroids(not(local_glia_locations),:);

    % Calculating the distances and converting from pixels to micrometers
    NtoN = pdist(neuronCentroids)*4620/12006;
    GtoG = pdist(gliaCentroids)*4620/12006;
    NtoC = (sqrt(neuronCentroids(:,1).^2 + ...
        neuronCentroids(:,2).^2))'*4620/12006;
    GtoC = (sqrt(gliaCentroids(:,1).^2 + ...
        gliaCentroids(:,2).^2))'*4620/12006;
    NtoNC = [reshape( ...
        pdist2(central_neurons, non_central_neurons).',1,[]), ...
        pdist(central_neurons)]*4620/12006;
    GtoGC = [reshape( ...
        pdist2(central_glia, non_central_glia).',1,[]), ...
        pdist(central_glia)]*4620/12006;
    NtoNL = [reshape( ...
        pdist2(local_neurons, non_local_neurons).',1,[]), ...
        pdist(local_neurons)]*4620/12006;
    GtoGL = [reshape( ...
        pdist2(local_glia, non_local_glia).',1,[]), ...
        pdist(local_glia)]*4620/12006;
    NtoGL = pdist2(neuronCentroids, gliaCentroids)*4620/12006;
    GtoNL = pdist2(gliaCentroids, neuronCentroids)*4620/12006;

    % Removing distance to itself (==0)
    NtoN(NtoN == 0) = [];
    GtoG(GtoG == 0) = [];
    NtoC(NtoC == 0) = [];
    GtoC(GtoC == 0) = [];
    NtoNC(NtoNC == 0) = [];
    GtoGC(GtoGC == 0) = [];
    NtoNL(NtoNL == 0) = [];
    GtoGL(GtoGL == 0) = [];
%     NtoGL(NtoGL == 0) = [];
%     GtoNL(GtoNL == 0) = [];

    % Do not need to save distances that are too long
    NtoNC(NtoNC >= 4620/4) = [];
    GtoGC(GtoGC >= 4620/4) = [];
    NtoNL(NtoNL >= 100) = [];
    GtoGL(GtoGL >= 100) = [];
    NtoGL(NtoGL >= 100) = 0;
    GtoNL(GtoNL >= 100) = 0;

end

