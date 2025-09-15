
% ----------------- Parameter Initialization -----------------
num_neurons = 550;  % Number of neurons

radius_well_mm = 7.75;  % Radius of the well in mm
radius_well_um = radius_well_mm * 1000;  % Convert to micrometers (µm)
cell_diameter_um = 14;  % Cell diameter in micrometers (µm)
critical_distance_um = 30;  % Increased Critical distance (e.g., 50 µm instead of 26 µm)

% Define conversion factor from µm to px
px_per_um = 1104 / 878;  % Conversion factor from µm to pixels

% Define the neuron radius in pixels
neuron_radius_um = cell_diameter_um / 2;
neuron_radius_px = neuron_radius_um * px_per_um;  % Convert to pixels

% ----------------- Generate Random Neuron Positions -----------------
neurons_px = [];  % Initialize the array of neuron positions

while length(neurons_px) < num_neurons
    % Randomly generate neuron positions
    theta = 2 * pi * rand(1, 1);
    r = radius_well_um * sqrt(rand(1, 1));  % Uniformly distributed within the circle
    x_um = r * cos(theta);  % x-coordinate in µm
    y_um = r * sin(theta);  % y-coordinate in µm
    
    % Convert coordinates from µm to pixels
    x_px = x_um * px_per_um;  % Convert x to pixels
    y_px = y_um * px_per_um;  % Convert y to pixels

    % Check for overlap with existing neurons (only if neurons already exist)
    if isempty(neurons_px)  % No neurons yet, so skip overlap check
        neurons_px = [neurons_px; x_px, y_px];
    else
        distances = sqrt((neurons_px(:, 1) - x_px).^2 + (neurons_px(:, 2) - y_px).^2);  % Euclidean distance to all other neurons
        if all(distances >= neuron_radius_px)  % Check if all distances are greater than the radius (no overlap)
            % If no overlap, add this neuron to the list
            neurons_px = [neurons_px; x_px, y_px];
        end
    end
end


% Now neurons_px contains positions of non-overlapping neurons

% ----------------- Define Centroids for Distance Analysis -----------------
neuronCentroids_px = neurons_px;  % Neuron positions in pixels
gliaCentroids_px = zeros(0, 2);  % Empty gliaCentroids but with correct size

% ----------------- Call the Distance Analysis Function -----------------
[NtoN, GtoG, NtoC, GtoC, NtoNC, GtoGC, NtoNL, ~, ~, ~] = Distance_analysis(neuronCentroids_px, gliaCentroids_px);
% Module 3 - Cluster Analysis Neurons


% Run the function, but collect graph data separately
[clusterData_neurons] = cluster_analysis_24well(neuronCentroids_px, critical_distance_um);

% % % ----------------- Identify Non-Local Neurons -----------------
% % well_radius_px = 12006 / 2;  % Given image size (radius in pixels)
% % edge_threshold_px = 100 * 12006 / 4620;  % Convert 100 µm to pixels
% % non_local_neurons = sqrt(neuronCentroids_px(:, 1).^2 + neuronCentroids_px(:, 2).^2) >= (well_radius_px - edge_threshold_px);
% % local_neurons = ~non_local_neurons;
% % 
% % % Extract coordinates of local and non-local neurons
% % local_x_px = neuronCentroids_px(local_neurons, 1);
% % local_y_px = neuronCentroids_px(local_neurons, 2);
% % non_local_x_px = neuronCentroids_px(non_local_neurons, 1);
% % non_local_y_px = neuronCentroids_px(non_local_neurons, 2);
%% NtoNL

% ----------------- Visualization of NtoNL Distances - Local Neuron Density -----------------
norm = input("Normalize, Yes or No? \n", "s");
if isempty(norm)
    norm = 'No';
end
n = '';  % Default naming suffix

% Define bin width
binwidth = 5; % Bin width in µm

% Load the NtoNL distances from the random neuron distribution
NtoNL = NtoNL(NtoNL < 500); % Keep only distances < 100 µm
NtoNL(NtoNL == 0) = []; % Remove zero distances (self-distance)

% Plot histogram
figure()
h = histogram(NtoNL, 'BinWidth', binwidth, 'FaceColor', 'g', 'EdgeColor', 'k');
hold on

% Normalize based on area (approximation)
h.BinCounts = h.BinCounts ./ (2 * pi * (h.BinEdges(1:end-1) + h.BinWidth) * h.BinWidth);

% Normalize if selected
if norm == "Yes"
    h.Normalization = 'probability';
    n = 'Normalized';
end

% Normalize based on number of neurons
h.BinCounts = h.BinCounts / length(NtoNL);

% Titles and labels
if norm == "Yes"
    title(['Random Neuron Distribution: Local NtoNL Distances (Normalized) - ', num2str(num_neurons), ' Neurons']);
else
    title(['Random Neuron Distribution: Local NtoNL Distances - ', num2str(num_neurons), ' Neurons']);
end
xlabel('Distance (\mu{\itm})');
ylabel('Local neuronal density');
legend("NtoNL Distances", 'Location', 'best');

% Set fixed y-axis maximum to 0.025
ylim([0, 0.04]);

% Define savepath for figure saving
savepath = 'E:\Angelica\Results\Random';  % Change this to your desired folder path

% Save the figure in the defined folder path
saveas(gcf, fullfile(savepath, strcat('Random_NtoNL_Distance', n, '.png')))
saveas(gcf, fullfile(savepath, strcat('Random_NtoNL_Distance', n, '.fig')))


% ----------------- Save Histogram Data to Excel -----------------
h_normalized = h.BinCounts / sum(h.BinCounts);
% Save normalized histograms to Excel file
if norm == "Yes"
% Dataset 1
    dataToSave = [h.BinEdges(1:end-1)', h_normalized'];
    columnTitles = {'Distance_um', 'NeuronDensity'};
% Construct the full path for the Excel file in the same directory as the figure
    excelFilePath = fullfile(savepath, 'Random_NtoNL_Distance.xlsx');
% Write data and column titles to Excel file on separate sheets
    writematrix(dataToSave, excelFilePath, 'Sheet', 'Dataset1', 'Range', 'A2');
    writecell(columnTitles, excelFilePath, 'Sheet', 'Dataset1', 'Range', 'A1');


end



%% 

%% NtoNC

% ----------------- Visualization of NtoNL Distances - Local Neuron Density -----------------
norm = input("Normalize, Yes or No? \n", "s");
if isempty(norm)
    norm = 'No';
end
n = '';  % Default naming suffix

% Define bin width
binwidth = 20; % Bin width in µm

% Load the NtoNL distances from the random neuron distribution
NtoNC = NtoNC(NtoNC < 4620/4 - mod(4620/4,binwidth)); 
NtoNC(NtoNC == 0) = []; % Remove zero distances (self-distance)

% Plot histogram
figure()
h = histogram(NtoNC, 'BinWidth', binwidth, 'FaceColor', 'g', 'EdgeColor', 'k');
hold on

% Normalize based on area (approximation)
h.BinCounts = h.BinCounts ./ (2 * pi * (h.BinEdges(1:end-1) + h.BinWidth) * h.BinWidth);

% Normalize if selected
if norm == "Yes"
    h.Normalization = 'probability';
    n = 'Normalized';
end

% Normalize based on number of neurons
h.BinCounts = h.BinCounts / length(NtoNC);

% Titles and labels
if norm == "Yes"
    title(['Random Neuron Distribution: Local NtoNC Distances (Normalized) - ', num2str(num_neurons), ' Neurons']);
else
    title(['Random Neuron Distribution: Local NtoNC Distances - ', num2str(num_neurons), ' Neurons']);
end
xlabel('Distance (\mu{\itm})');
ylabel('Local neuronal density');
legend("NtoNL Distances", 'Location', 'best');

% Set fixed y-axis maximum to 0.025
ylim([0, 0.03]);

% Define savepath for figure saving
savepath = 'E:\Angelica\Results\Random';  % Change this to your desired folder path

% Save the figure in the defined folder path
saveas(gcf, fullfile(savepath, strcat('Random_NtoNC_Distance', n, '.png')))
saveas(gcf, fullfile(savepath, strcat('Random_NtoNC_Distance', n, '.fig')))


% ----------------- Save Histogram Data to Excel -----------------
h_normalized = h.BinCounts / sum(h.BinCounts);
% Save normalized histograms to Excel file
if norm == "Yes"
% Dataset 1
    dataToSave = [h.BinEdges(1:end-1)', h_normalized'];
    columnTitles = {'Distance_um', 'NeuronDensity'};
% Construct the full path for the Excel file in the same directory as the figure
    excelFilePath = fullfile(savepath, 'Random_NtoNC_Distance.xlsx');
% Write data and column titles to Excel file on separate sheets
    writematrix(dataToSave, excelFilePath, 'Sheet', 'Dataset1', 'Range', 'A2');
    writecell(columnTitles, excelFilePath, 'Sheet', 'Dataset1', 'Range', 'A1');


end


%% 

%% NtoN

% ----------------- Visualization of NtoN Distances - Local Neuron Density -----------------
norm = input("Normalize, Yes or No? \n", "s");
if isempty(norm)
    norm = 'No';
end
n = '';  % Default naming suffix

% Define bin width
binwidth = 50; % Bin width in µm

% Load the NtoN distances from the random neuron distribution

NtoN (NtoN == 0) = []; % Remove zero distances (self-distance)

% Plot histogram
figure()
h = histogram(NtoN, 'BinWidth', binwidth, 'FaceColor', 'g', 'EdgeColor', 'k');
hold on

% Normalize based on area (approximation)
h.BinCounts = h.BinCounts ./ (2 * pi * (h.BinEdges(1:end-1) + h.BinWidth) * h.BinWidth);

% Normalize if selected
if norm == "Yes"
    h.Normalization = 'probability';
    n = 'Normalized';
end

% Normalize based on number of neurons
h.BinCounts = h.BinCounts / length(NtoN);

% Titles and labels
if norm == "Yes"
    title(['Random Neuron Distribution: Local NtoN Distances (Normalized) - ', num2str(num_neurons), ' Neurons']);
else
    title(['Random Neuron Distribution: Local NtoN Distances - ', num2str(num_neurons), ' Neurons']);
end
xlabel('Distance (\mu{\itm})');
ylabel('Local neuronal density');
legend("NtoNL Distances", 'Location', 'best');

% Set fixed y-axis maximum to 0.025
%ylim([0, 0.03]);

% Define savepath for figure saving
savepath = 'E:\Angelica\Results\Random';  % Change this to your desired folder path

% Save the figure in the defined folder path
saveas(gcf, fullfile(savepath, strcat('Random_NtoN_Distance', n, '.png')))
saveas(gcf, fullfile(savepath, strcat('Random_NtoN_Distance', n, '.fig')))


% ----------------- Save Histogram Data to Excel -----------------
h_normalized = h.BinCounts / sum(h.BinCounts);
% Save normalized histograms to Excel file
if norm == "Yes"
% Dataset 1
    dataToSave = [h.BinEdges(1:end-1)', h_normalized'];
    columnTitles = {'Distance_um', 'NeuronDensity'};
% Construct the full path for the Excel file in the same directory as the figure
    excelFilePath = fullfile(savepath, 'Random_NtoN_Distance.xlsx');
% Write data and column titles to Excel file on separate sheets
    writematrix(dataToSave, excelFilePath, 'Sheet', 'Dataset1', 'Range', 'A2');
    writecell(columnTitles, excelFilePath, 'Sheet', 'Dataset1', 'Range', 'A1');


end

%% 
%% NtoC

% ----------------- Visualization of NtoNL Distances - Local Neuron Density -----------------
norm = input("Normalize, Yes or No? \n", "s");
if isempty(norm)
    norm = 'No';
end
n = '';  % Default naming suffix

% Define bin width
binwidth = 50; % Bin width in µm

% Load the NtoNL distances from the random neuron distribution

NtoC(NtoC == 0) = []; % Remove zero distances (self-distance)

% Plot histogram
figure()
h = histogram(NtoC, 'BinWidth', binwidth, 'FaceColor', 'g', 'EdgeColor', 'k');
hold on

% Normalize based on area (approximation)
h.BinCounts = h.BinCounts ./ (2 * pi * (h.BinEdges(1:end-1) + h.BinWidth) * h.BinWidth);

% Normalize if selected
if norm == "Yes"
    h.Normalization = 'probability';
    n = 'Normalized';
end

% Normalize based on number of neurons
h.BinCounts = h.BinCounts / length(NtoC);

% Titles and labels
if norm == "Yes"
    title(['Random Neuron Distribution: NtoC Distances (Normalized) - ', num2str(num_neurons), ' Neurons']);
else
    title(['Random Neuron Distribution: NtoC Distances - ', num2str(num_neurons), ' Neurons']);
end
xlabel('Distance (\mu{\itm})');
ylabel('Local neuronal density');
legend("NtoNL Distances", 'Location', 'best');

% Set fixed y-axis maximum to 0.025
%ylim([0, 0.03]);

% Define savepath for figure saving
savepath = 'E:\Angelica\Results\Random';  % Change this to your desired folder path

% Save the figure in the defined folder path
saveas(gcf, fullfile(savepath, strcat('Random_NtoC_Distance', n, '.png')))
saveas(gcf, fullfile(savepath, strcat('Random_NtoC_Distance', n, '.fig')))


% ----------------- Save Histogram Data to Excel -----------------
h_normalized = h.BinCounts / sum(h.BinCounts);
% Save normalized histograms to Excel file
if norm == "Yes"
% Dataset 1
    dataToSave = [h.BinEdges(1:end-1)', h_normalized'];
    columnTitles = {'Distance_um', 'NeuronDensity'};
% Construct the full path for the Excel file in the same directory as the figure
    excelFilePath = fullfile(savepath, ['Random_NtoC_Distance.xlsx']);
% Write data and column titles to Excel file on separate sheets
    writematrix(dataToSave, excelFilePath, 'Sheet', 'Dataset1', 'Range', 'A2');
    writecell(columnTitles, excelFilePath, 'Sheet', 'Dataset1', 'Range', 'A1');


end
%% 




