%% Load settings
MinPeakDistancesce = 5;      % 5 default
MinPeakDistance = 3;          % 3 default
threshold_peak = 3;           % 3 default
synchronous_frames = 2;       % 2 default
sce_n_cells_threshold = 10;
percentile = NaN;
minithreshold=0.1;
SG_window = 9 ;

%% Load current data for motion correction and colorsubstraction
load(strcat(path, 'Fall.mat'), 'ops');
fileExists_colorcellnew = isfile(strcat(path ,'colorcellnew.mat'));
fileExists_colorcell = isfile(strcat(path ,'colorcell.mat'));
if fileExists_colorcellnew
    load (strcat(path ,'colorcellnew.mat'))
elseif fileExists_colorcell
    load (strcat(path ,'colorcell.mat'))
end
% load("E:\Data\Aurelie\data\nocues\444119\220919_plane0\colorcell_registration.mat");%colorcell from Solene;%for444119 220919plane0
 Tr1b = double(F);

if exist('colorcell','var')
    mask_to_keep = colorcell < 7;
    Tr1b = Tr1b(mask_to_keep, :);
    colorcell = colorcell(mask_to_keep);
end
%% ===== DÉTECTION DES BAD FRAMES =====

corrXY = ops.corrXY;

% Approche robuste : Déviation par rapport à la tendance locale
rolling_median = movmedian(corrXY, 300);
deviation = corrXY - rolling_median;

% Bad frames = celles qui dévient fortement vers le bas
sigma_dev = std(deviation(deviation < 0));
seuil_bad = -3 * sigma_dev;
bad_frames = deviation < seuil_bad;
bad_frames = conv(double(bad_frames), [1 1 1], 'same') > 0;% to be safe we remove previous and next frames

% fprintf('Bad frames détectées : %d (%.2f%%)\n', ...
%     sum(bad_frames), 100*sum(bad_frames)/length(corrXY));

% Bad frames sans mouvement
bad_frames_no_movement = bad_frames & (speed' < 2);
n_bad_no_move = sum(bad_frames_no_movement);
% fprintf('Bad frames avec speed < 2 cm/s : %d (%.2f%%)\n', ...
%     n_bad_no_move, 100*n_bad_no_move/length(corrXY));

% % Visualisation
% figure;
% subplot(2,1,1);
% plot(corrXY, 'b'); hold on;
% plot(rolling_median, 'r-', 'LineWidth', 2);
% plot(find(bad_frames), corrXY(bad_frames), 'rx', 'MarkerSize', 8, 'DisplayName', 'Bad Frames (tous)');
% plot(find(bad_frames_no_movement), corrXY(bad_frames_no_movement), 'ko', 'MarkerSize', 8, 'DisplayName', 'Bad Frames sans mouvement');
% title('Bad Frames (Détection par Déviation Locale)');
% legend;
% ylabel('ops.corrXY');
% 
% subplot(2,1,2);
% plot(speed, 'g'); hold on;
% yline(2, 'r--', 'Speed Threshold');
% scatter(find(bad_frames_no_movement), speed(bad_frames_no_movement), 50, 'ko', 'filled', 'DisplayName', 'Bad frames (speed < 2)');
% title('Vitesse de la Souris');
% ylabel('Speed (cm/s)');
% xlabel('Frames');
% legend;
% 
% fprintf('\n=== Résumé Bad Frames ===\n');
% fprintf('Total bad frames : %d\n', sum(bad_frames));
% fprintf('Bad frames avec mouvement (speed >= 2) : %d\n', sum(bad_frames & (speed' >= 2)));
%fprintf('Bad frames sans mouvement (speed < 2) : %d\n', n_bad_no_move);

% ===== PRÉPARATION DES DONNÉES =====

speedsm = smoothdata(speed, 'gaussian', 50);
[NCell, Nz] = size(Tr1b);
% fprintf('Ncells = %d\n', NCell);

% Bleaching correction
Tr1b_clean = Tr1b;
Tr1b_clean(:, bad_frames) = NaN;
Tr1b_clean = fillmissing(Tr1b_clean, 'linear', 2, 'EndValues', 'nearest'); %interpolation
Tr1b = sgolayfilt(Tr1b_clean', 3, SG_window)';

% ===== NORMALISATION dF/F AVEC FILTRAGE DES BAD FRAMES =====
window_size = 1000;
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
%%%% 1ere reduction de dimension 



% ===== CRÉER DES WINDOWS VALIDES (SANS BAD FRAMES) =====
WinRest = find((speedsm <= 1) & ~bad_frames' & (speedsm > 0));
% WinActive = find((speedsm > 2) & ~bad_frames');
WinActive = find(speedsm > 2);

% ===== DÉTECTION DES TRANSIENTS =====
Raster = zeros(NCell, Nz);
Acttmp2 = cell(1, NCell);
ampli = cell(1, NCell);
th = zeros(1, NCell);

for i = 1:NCell
    % Seuil calculé SANS les bad frames et SANS le mouvement
    mad_trace = mad(Tr1b(i, WinRest), 1);
    th(i) = threshold_peak * mad_trace * 1.4826;
    
    % Détecter les pics sur la trace complète
    [amplitude, locs] = findpeaks(Tr1b(i, :), 'MinPeakProminence', th(i), 'MinPeakDistance', MinPeakDistance);
    
    % Filtrer : enlever les pics qui tombent sur des bad frames
    valid_mask = ~bad_frames(locs)';
    locs = locs(valid_mask);
    amplitude = amplitude(valid_mask);
    
    % Filtrer les pics en mouvement
    valeurs_identiques = intersect(locs, WinActive);
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
% nombre_transients_par_cellule = cellfun(@length, Acttmp2);
% seuil_frequence = prctile(nombre_transients_par_cellule, 99);
% cellules_hyperactives_idx = find(nombre_transients_par_cellule > seuil_frequence);
% if ~isempty(cellules_hyperactives_idx)
%     Raster(cellules_hyperactives_idx, :) = 0;
% end

% ===== DÉTECTION DES SCE (Synchronous Events) =====
%%%%shuffling to find threshold for number of cell for sce detection
% MActsh = zeros(1,Nz-synchronous_frames);   
% Rastersh=zeros(NCell,Nz);   
% NShfl=100;
% Sumactsh=zeros(Nz-synchronous_frames,NShfl);   
% for n=1:NShfl
% 
%         for c=1:NCell
%             k = randi(Nz-length(WinActive));
%             Rastersh(c,:)= circshift(Raster(c,:),k,2);
%         end
% 
%         for i=1:Nz-synchronous_frames   %need to use WinRest???
%             MActsh(i) = sum(max(Rastersh(:,i:i+synchronous_frames),[],2));
%         end
% 
%     Sumactsh(:,n)=MActsh;
% end
% 
% percentile = 95; % Calculate the 5% highest point or 99
% sce_n_cells_threshold = prctile(Sumactsh, percentile,"all");


SumAct = sum(Raster, 1);
[pks, locs] = findpeaks(SumAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);

% Filtrer les SCE : enlever celles qui tombent sur des bad frames
sces_valid_mask = ~bad_frames(locs)';
locs = locs(sces_valid_mask);
pks = pks(sces_valid_mask);

% Filtrer les SCE aberrantes (outliers)
% absolute_threshold = 100;
% TF_relative = isoutlier(pks, "percentiles", [0 99]);
% TF_absolute = pks > absolute_threshold;
% TF_combined = TF_relative & TF_absolute;
% idx_to_remove = find(TF_combined);
% if ~isempty(idx_to_remove)
%     Raster(:, locs(idx_to_remove)) = 0;
% end

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

% Paramètres
nShuffles = 100;  % 50 est généralement suffisant pour stabiliser le seuil
alpha = 95;      % Seuil de confiance (95ème percentile)

% 1. PCA sur données réelles
warning('off', 'all');
[~, score, latent_real] = pca(Race'); 
nPCs = length(latent_real);
[nTime, nNeurons] = size(Race');

% Matrice pour stocker les valeurs propres de chaque shuffle
null_latents = zeros(nPCs, nShuffles);

% 2. Boucle de Shuffles (peut prendre quelques secondes/minutes)
% fprintf('Lancement de %d shuffles...\n', nShuffles);

for s = 1:nShuffles
    % Création d'une matrice aléatoire (Circular Shift)
    Race_Shuffled = zeros(nTime, nNeurons);
    for i = 1:nNeurons
        shift = randi(nTime); 
        % Décalage circulaire pour garder l'autocorrélation temporelle de la cellule
        Race_Shuffled(:, i) = circshift(Race(i, :)', shift);
    end
    
    % PCA sur le bruit
    [~, ~, latent_null] = pca(Race_Shuffled);
    null_latents(:, s) = latent_null;
end
warning('on', 'all');

% 3. Calcul du seuil statistique (Parallel Analysis)
% Pour chaque PC, on regarde quelle est la valeur propre max du bruit dans 95% des cas
threshold_curve = prctile(null_latents, alpha, 2); 

% 4. Comparaison : Réel vs Bruit
% On garde tant que la vraie valeur propre est supérieure au seuil de bruit
keep_mask = latent_real > threshold_curve;

% On cherche le PREMIER moment où le signal réel passe SOUS le bruit
first_crossing = find(latent_real < threshold_curve, 1, 'first');

if isempty(first_crossing)
    % Si ça ne croise jamais (cas impossible mathématiquement sauf bug)
    nPC_Final = nPCs; 
else
    % On prend juste avant le croisement
    % nPC_Final = first_crossing - 1;
    nPC_Final = first_crossing + 5;
end

%fprintf('Nombre optimal corrigé : %d\n', nPC_Final);


% Visualisation pour vérifier (Optionnel mais recommandé)
% figure;
% plot(1:50, latent_real(1:50), '-o', 'LineWidth', 2, 'DisplayName', 'Données Réelles'); hold on;
% plot(1:50, threshold_curve(1:50), '--r', 'LineWidth', 2, 'DisplayName', 'Seuil Bruit (95%)');
% xline(nPC_Final, 'k:', ['Cutoff: ' num2str(nPC_Final)]);
% legend;
% xlabel('Composante Principale'); ylabel('Variance (Eigenvalue)');
% title('Parallel Analysis avec Shuffles');
% grid on;
% 
% fprintf('Nombre optimal de PC (Parallel Analysis) : %d\n', nPC_Final);
nPC_Final = min (nPC_Final,NCell/10);

% Paramètres
corr_threshold = 0.70; % Si corr > 0.7, on considère que c'est un artefact global
max_start_pc = 4;      % On ne rejette pas au-delà de la PC4 par sécurité

% 2. Calcul de la Trace Moyenne de la population (Signal Global)
% Moyenne de l'activité de toutes les cellules à chaque instant T
global_signal = mean(Race, 1)'; % (nTime x 1)

% 3. Boucle de décision automatique
start_PC = 1; % Par défaut on commence à 1

for i = 1:max_start_pc
    % Calcul de la corrélation absolue entre la PC(i) et le signal global
    r = corr(score(:, i), global_signal);
    
    % fprintf('PC %d vs Global Signal : Correlation = %.2f ', i, r);
    
    if abs(r) > corr_threshold
        % fprintf('-> REJETÉE (Artefact probable)\n');
        start_PC = i + 1; % On décale le début
    else
        % fprintf('-> GARDÉE (Signal spécifique)\n');
        break; % Dès qu'une PC est "propre", on s'arrête et on garde tout le reste
    end
end
% % 
start_PC=2;
nPC_Final=100;

Race_For_Clustering = score(:, start_PC:nPC_Final)';
