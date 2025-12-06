%% Load settings
MinPeakDistancesce = 5;      % 5 default
MinPeakDistance = 3;          % 3 default
threshold_peak = 3;           % 3 default
synchronous_frames = 2;       % 2 default
sce_n_cells_threshold = 5;
percentile = NaN;

%% Load current data
load(strcat(path, 'Fall.mat'), 'ops');

% ===== DÉTECTION DES BAD FRAMES =====
corrXY = ops.corrXY;

% Approche robuste : Déviation par rapport à la tendance locale
rolling_median = movmedian(corrXY, 300);
deviation = corrXY - rolling_median;

% Bad frames = celles qui dévient fortement vers le bas
sigma_dev = std(deviation(deviation < 0));
seuil_bad = -3 * sigma_dev;
bad_frames = deviation < seuil_bad;

fprintf('Bad frames détectées : %d (%.2f%%)\n', ...
    sum(bad_frames), 100*sum(bad_frames)/length(corrXY));

% Bad frames sans mouvement
bad_frames_no_movement = bad_frames & (speed' < 2);
n_bad_no_move = sum(bad_frames_no_movement);
fprintf('Bad frames avec speed < 2 cm/s : %d (%.2f%%)\n', ...
    n_bad_no_move, 100*n_bad_no_move/length(corrXY));

% Visualisation
figure;
subplot(2,1,1);
plot(corrXY, 'b'); hold on;
plot(rolling_median, 'r-', 'LineWidth', 2);
plot(find(bad_frames), corrXY(bad_frames), 'rx', 'MarkerSize', 8, 'DisplayName', 'Bad Frames (tous)');
plot(find(bad_frames_no_movement), corrXY(bad_frames_no_movement), 'ko', 'MarkerSize', 8, 'DisplayName', 'Bad Frames sans mouvement');
title('Bad Frames (Détection par Déviation Locale)');
legend;
ylabel('ops.corrXY');

subplot(2,1,2);
plot(speed, 'g'); hold on;
yline(2, 'r--', 'Speed Threshold');
scatter(find(bad_frames_no_movement), speed(bad_frames_no_movement), 50, 'ko', 'filled', 'DisplayName', 'Bad frames (speed < 2)');
title('Vitesse de la Souris');
ylabel('Speed (cm/s)');
xlabel('Frames');
legend;

fprintf('\n=== Résumé Bad Frames ===\n');
fprintf('Total bad frames : %d\n', sum(bad_frames));
fprintf('Bad frames avec mouvement (speed >= 2) : %d\n', sum(bad_frames & (speed' >= 2)));
fprintf('Bad frames sans mouvement (speed < 2) : %d\n', n_bad_no_move);

% ===== PRÉPARATION DES DONNÉES =====
Tr1b = double(F);
speedsm = smoothdata(speed, 'gaussian', 50);
[NCell, Nz] = size(Tr1b);
fprintf('Ncells = %d\n', NCell);

% Bleaching correction
Tr1b = sgolayfilt(Tr1b', 3, 9)';

% ===== NORMALISATION dF/F AVEC FILTRAGE DES BAD FRAMES =====
window_size = 500;
percentile_value = 5;
num_blocks = ceil(Nz / window_size);

for n = 1:NCell
    trace = Tr1b(n, :);
    
    % Masquer TEMPORAIREMENT les bad frames pour le calcul de F0
    trace_masked = trace;
    trace_masked(bad_frames') = nan;
    
    F0 = nan(Nz, 1);
    % Calcul du percentile bloc par bloc
    for i = 1:num_blocks
        start_idx = (i-1) * window_size + 1;
        end_idx = min(i * window_size, Nz);
        F0(start_idx:end_idx) = prctile(trace_masked(start_idx:end_idx), percentile_value);
    end
    
    F0 = movmedian(F0, window_size, 'omitnan');
    F0 = smoothdata(F0, 1, "gaussian", window_size/2);
    
    % dF/F sur la trace ORIGINALE (pas de nan)
    Tr1b(n, :) = (trace - F0') ./ F0';
end

% ===== CRÉER DES WINDOWS VALIDES (SANS BAD FRAMES) =====
WinRest_valid = find((speedsm <= 2) & ~bad_frames');
WinActive_valid = find((speedsm > 2) & ~bad_frames');

% ===== DÉTECTION DES TRANSIENTS =====
Raster = zeros(NCell, Nz);
Acttmp2 = cell(1, NCell);
ampli = cell(1, NCell);
th = zeros(1, NCell);

for i = 1:NCell
    % Seuil calculé SANS les bad frames et SANS le mouvement
    mad_trace = mad(Tr1b(i, WinRest_valid), 1);
    th(i) = threshold_peak * mad_trace * 1.4826;
    
    % Détecter les pics sur la trace complète
    [amplitude, locs] = findpeaks(Tr1b(i, :), 'MinPeakProminence', th(i), 'MinPeakDistance', MinPeakDistance);
    
    % Filtrer : enlever les pics qui tombent sur des bad frames
    valid_mask = ~bad_frames(locs)';
    locs = locs(valid_mask);
    amplitude = amplitude(valid_mask);
    
    % Filtrer les pics en mouvement
    valeurs_identiques = intersect(locs, WinActive_valid);
    [locs_sans_ide, idx] = setdiff(locs(:), valeurs_identiques);
    ampli_sans_ide = amplitude(idx);
    Acttmp2{i} = locs_sans_ide;
end

% Créer Raster avec les transients valides
for i = 1:NCell
    if ~isempty(Acttmp2{i})
        Raster(i, Acttmp2{i}) = 1;
    end
end

% ===== ÉLIMINER LES CELLULES HYPERACTIVES =====
nombre_transients_par_cellule = cellfun(@length, Acttmp2);
seuil_frequence = prctile(nombre_transients_par_cellule, 98);
cellules_hyperactives_idx = find(nombre_transients_par_cellule > seuil_frequence);
if ~isempty(cellules_hyperactives_idx)
    Raster(cellules_hyperactives_idx, :) = 0;
end

% ===== DÉTECTION DES SCE (Synchronous Events) =====
SumAct = sum(Raster, 1);
[pks, locs] = findpeaks(SumAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);

% Filtrer les SCE : enlever celles qui tombent sur des bad frames
sces_valid_mask = ~bad_frames(locs)';
locs = locs(sces_valid_mask);
pks = pks(sces_valid_mask);

% Filtrer les SCE aberrantes (outliers)
absolute_threshold = 100;
TF_relative = isoutlier(pks, "percentiles", [0 95]);
TF_absolute = pks > absolute_threshold;
TF_combined = TF_relative & TF_absolute;
idx_to_remove = find(TF_combined);
if ~isempty(idx_to_remove)
    Raster(:, locs(idx_to_remove)) = 0;
end

% Recalculer après filtrage
SumAct = sum(Raster, 1);
[pks, TRace] = findpeaks(SumAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);

NRace = length(TRace);
fprintf('nSCE après filtrage bad_frames : %d\n', NRace);

% ===== CRÉER RACE MATRIX (Participation des cellules aux SCE) =====
MAct = zeros(1, Nz - synchronous_frames);
for i = 1:Nz - synchronous_frames
    MAct(i) = sum(max(Raster(:, i:i+synchronous_frames), [], 2));
end

fprintf('Sum transient: %d\n', sum(MAct));

Race = false(NCell, NRace);
RasterRace = zeros(NCell, Nz);
for i = 1:NRace
    Race(:, i) = max(Raster(:, TRace(i)-1:TRace(i)+2), [], 2);
    RasterRace(Race(:, i) == 1, TRace(i)) = 1;
end

Race = double(Race);

% ===== PCA =====
warning('off', 'all');
[~, score, ~, ~, explained] = pca(Race');
warning('on', 'all');

cumVar = cumsum(explained);
nPC = find(cumVar >= 80, 1);
Race_For_Clustering = score(:, 3:30)';

fprintf('\nPCA : %d PCs pour 80%% de variance\n', nPC);
