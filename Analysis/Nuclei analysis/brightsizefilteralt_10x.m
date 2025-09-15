function [Opt_LBW, Data, threshold] = brightsizefilteralt_10x(LBW_Neurons, Dapi)
    %% Initialize output and calculate connected components
    Opt_LBW = LBW_Neurons;  % Copy input to output (binary image of neurons)
    CC = bwconncomp(Opt_LBW);  % Find connected components in binary image
    
    % Calculate region properties: Area, Mean Intensity, Pixel Values (intensities)
    stats = regionprops('table', CC, Dapi, 'Area', 'MeanIntensity', 'PixelValues');
    
    % Calculate standard deviation of pixel intensities for each region
    sigma = zeros(length(stats.PixelValues), 1);
    for i = 1:length(stats.PixelValues)
        sigma(i) = std(double(stats.PixelValues{i}));
    end
    
    %% Clustering based on region area and intensity variation (standard deviation)
    % Scale area by a factor of 4 for 10x magnification (adjust from 20x to 10x)
    scaled_area = stats.Area / 4;
    
    % Perform k-medoids clustering using the scaled area and intensity variation
    [idx, C] = kmedoids([scaled_area sigma], 2, 'Distance', 'seuclidean', 'Replicates', 5);
    
    % Determine which cluster corresponds to higher intensity variation
    if C(1, 2) > C(2, 2)
        Idx_intensity = (idx == 1);  % Cluster 1 has higher variation in intensity
    else
        Idx_intensity = (idx == 2);  % Cluster 2 has higher variation in intensity
    end
    
    %% Remove regions with higher intensity variation (likely noise)
    % Set pixels corresponding to these regions to false (remove them)
    Opt_LBW(cat(1, CC.PixelIdxList{Idx_intensity})) = false;
    
    %% Find local minima in the area distribution to set a size threshold
    % Create a histogram of region areas (50 bins)
    [N, X] = hist(stats.Area, 50);
    
    % Find local minima in the histogram
    [Minima, P] = islocalmin(N);
    
    % Set the threshold as the area corresponding to the maximum local minimum
    threshold = X(P == max(P));
    threshold = threshold(1);  % In case there are multiple, take the first one
    
    % Cap the threshold at 550 (scaled for 20x, but we adjust this for 10x)
    if threshold > 550 / 4
        threshold = 550 / 4;  % Adjust the threshold for 10x magnification
    end
    
    %% Filter out small regions based on the area threshold
    % Identify regions whose area is below the threshold
    Idx_area = (stats.Area < threshold);
    
    % Remove regions that are too small
    Opt_LBW(cat(1, CC.PixelIdxList{Idx_area})) = false;
    
    %% Output Data
    % Store region area, intensity variation, and cluster index for further analysis
    Data = [stats.Area sigma idx];
    
end