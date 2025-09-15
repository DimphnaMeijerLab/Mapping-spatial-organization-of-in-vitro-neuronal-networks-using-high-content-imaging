%% AIS plotting & export (Div_7_segmented / B08)
clc; close all; clear;

%% ===== Settings =====
rootPath = 'E:\Cellpose\HP\Div_7_segmented\B08';   % adjust if needed
dataset  = 'new_conversion_correct';
well     = 'B08';
savePath  = uigetdir(pwd, 'Select output folder');
if savePath == 0, error('No output folder selected'); end

% Pixel-to-micron conversion
pixelSize = 1;   % ≈ 0.401 µm/px

%% ===== Load data =====
matFile = fullfile(rootPath, dataset, well, 'AISResults.mat');
if ~isfile(matFile)
    error('Missing file: %s', matFile);
end

S = load(matFile);
if isfield(S,'AISResults')
    Ais = S.AISResults;
elseif isfield(S,'AIS_Results')
    Ais = S.AIS_Results;
else
    error('No AISResults struct found in %s', matFile);
end

%% ===== Convert values =====
% AIS lengths
AISlen = {};
if isfield(Ais,'AISlength')
    AISlen = cellfun(@(x) x * pixelSize, Ais.AISlength, 'UniformOutput', false);
end
rawLengths = vertcat(AISlen{:});

% Distances
dist = [];
if isfield(Ais,'Distances')
    dist = cellfun(@(x) x * pixelSize, Ais.Distances, 'UniformOutput', false);
    dist = vertcat(dist{:});
end

% Angles
ang = [];
if isfield(Ais,'Allangles') && ~isempty(Ais.Allangles) && isnumeric(Ais.Allangles)
    ang = toAnglesDegrees(Ais.Allangles(:,1));
elseif isfield(Ais,'Angles')
    ang = toAnglesDegrees(Ais.Angles);
end

% AIS count per neuron
aisCountPerNeuron = cellfun(@numel, Ais.AISlength);

%% ===== Figures =====
% 1. AIS length histogram
if ~isempty(rawLengths)
    f1 = figure('Color','w');
    histogram(rawLengths, 0:1:100);
    xlabel('AIS length (µm)');
    ylabel('Neuron count');
    title('AIS length distribution');
    set(gca,'FontSize',14,'LineWidth',1);
    saveas(f1, fullfile(savePath,'AIS_length_histogram.png'));
    close(f1);
end

% 2. AIS per neuron histogram
if ~isempty(aisCountPerNeuron)
    f2 = figure('Color','w');
    histogram(aisCountPerNeuron, 'BinMethod','integers');
    xlabel('AIS per neuron');
    ylabel('Neuron count');
    title('AIS per neuron distribution');
    set(gca,'FontSize',14,'LineWidth',1);
    saveas(f2, fullfile(savePath,'AIS_perNeuron_histogram.png'));
    close(f2);
end

% 3. Distance histogram
if ~isempty(dist)
    f3 = figure('Color','w');
    histogram(dist, 0:5:500);
    xlabel('Distance to soma (µm)');
    ylabel('AIS count');
    title('AIS distance distribution');
    set(gca,'FontSize',14,'LineWidth',1);
    saveas(f3, fullfile(savePath,'AIS_distance_histogram.png'));
    close(f3);
end

% 4. Angle histogram
if ~isempty(ang)
    f4 = figure('Color','w');
    histogram(ang, 0:10:360);
    xlabel('Angle (deg)');
    ylabel('AIS count');
    title('AIS angle distribution');
    set(gca,'FontSize',14,'LineWidth',1);
    saveas(f4, fullfile(savePath,'AIS_angle_histogram.png'));
    close(f4);

    % Polar histogram
    f5 = figure('Color','w');
    polarhistogram(deg2rad(ang), 18);
    title('AIS angle polar distribution');
    saveas(f5, fullfile(savePath,'AIS_angle_polar.png'));
    close(f5);
end

%% ===== Excel Export =====
outXLS = fullfile(savePath,'AIS_B08_Div7_results.xlsx');

% Raw lengths
if ~isempty(rawLengths)
    Tlen = table(rawLengths, 'VariableNames', {'AIS_length_um'});
    writetable(Tlen, outXLS, 'Sheet','RawLengths');
end

% AIS per neuron
if ~isempty(aisCountPerNeuron)
    vals = aisCountPerNeuron(:);
    edges = -0.5:1:(max(vals)+0.5);
    counts = histcounts(vals, edges)';
    mids   = edges(1:end-1)' + 0.5;
    perc   = counts / sum(counts) * 100;
    Tcount = table(mids, counts, perc, ...
        'VariableNames', {'AIS_per_neuron','Count','Percent'});
    writetable(Tcount, outXLS, 'Sheet','AISperNeuron');
end

% Distances
if ~isempty(dist)
    Tdist = table(dist, 'VariableNames', {'AIS_distance_um'});
    writetable(Tdist, outXLS, 'Sheet','Distances');
end

% Angles
if ~isempty(ang)
    Tang = table(ang, 'VariableNames', {'AIS_angle_deg'});
    writetable(Tang, outXLS, 'Sheet','Angles');

    edges = 0:10:360;
    counts = histcounts(ang, edges)';
    mids   = edges(1:end-1)' + 5;
    perc   = counts / sum(counts) * 100;
    TangHist = table(mids, counts, perc, ...
        'VariableNames', {'Angle_bin_mid_deg','Count','Percent'});
    writetable(TangHist, outXLS, 'Sheet','Angles_hist');
end

fprintf('✅ Saved Excel and PNG figures in: %s\n', savePath);

%% ===================== Helpers =====================
function aDeg = toAnglesDegrees(x)
    if isempty(x), aDeg = []; return; end
    a = double(x(:));
    if max(abs(a),[],'omitnan') <= pi*1.1
        a = rad2deg(a);
    end
    aDeg = mod(a, 360);
end
