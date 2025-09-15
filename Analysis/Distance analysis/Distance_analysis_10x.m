function [NtoN, NtoC, NtoNC, GtoGC, NtoNL, GtoGL] = Distance_analysis_10x(neuronCentroids, gliaCentroids)
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
        neuronCentroids(:,2).^2) < (19136/4);
    central_glia_locations = sqrt(gliaCentroids(:,1).^2 + ...
        gliaCentroids(:,2).^2) < (19136/4);
    central_neurons = neuronCentroids(central_neurons_locations,:);
    central_glia = gliaCentroids(central_glia_locations,:);
    non_central_neurons = neuronCentroids(not(central_neurons_locations),:);
    non_central_glia = gliaCentroids(not(central_glia_locations),:);

    % Find the cells with their centroid not closer than 100 um to the edge
    local_neurons_locations = sqrt(neuronCentroids(:,1).^2 + ...
        neuronCentroids(:,2).^2) < (19136/2-(100*1104/878));
    local_glia_locations = sqrt(gliaCentroids(:,1).^2 + ...
        gliaCentroids(:,2).^2) < (19136/2-(100*1104/878));
    local_neurons = neuronCentroids(local_neurons_locations,:);
    local_glia = gliaCentroids(local_glia_locations,:);
    non_local_neurons = neuronCentroids(not(local_neurons_locations),:);
    non_local_glia = gliaCentroids(not(local_glia_locations),:);

    % Calculating the distances and converting from pixels to micrometers
    NtoN = pdist(neuronCentroids)*878/1104;
    NtoC = (sqrt(neuronCentroids(:,1).^2 + ...
        neuronCentroids(:,2).^2))'*878/1104;
    NtoNC = [reshape( ...
        pdist2(central_neurons, non_central_neurons).',1,[]), ...
        pdist(central_neurons)]*878/1104;
    GtoGC = [reshape( ...
        pdist2(central_glia, non_central_glia).',1,[]), ...
        pdist(central_glia)]*878/1104;
    NtoNL = [reshape( ...
        pdist2(local_neurons, non_local_neurons).',1,[]), ...
        pdist(local_neurons)]*878/1104;
    GtoGL = [reshape( ...
        pdist2(local_glia, non_local_glia).',1,[]), ...
        pdist(local_glia)]*878/1104;

    % Removing distance to itself (==0)
    NtoN(NtoN == 0) = [];
    NtoC(NtoC == 0) = [];
    NtoNC(NtoNC == 0) = [];
    GtoGC(GtoGC == 0) = [];
    NtoNL(NtoNL == 0) = [];
    GtoGL(GtoGL == 0) = [];


    % Do not need to save distances that are too long
    NtoNC(NtoNC >= 15219/4) = [];
    GtoGC(GtoGC >= 15219/4) = [];
    NtoNL(NtoNL >= 100) = [];
    GtoGL(GtoGL >= 100) = [];


end

