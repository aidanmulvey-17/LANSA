%% LANSA Analysis Pipeline
% Principal component analysis (PCA) and partial least squares regression
% (PLSR) on the LANSA (Laminar Anatomy and Neural Systems Atlas) dataset.
%
% Mulvey et al., "Conserved Laminar-Cellular Architecture Predicts the
% Spectrolaminar Motif".
%
% Inputs
%   LANSA_processed_cyto.mat         Cytoarchitecture (included): laminar density of
%                            PV, CB, CR and NG cells in 11 areas, in the
%                            standardized laminar space (21 bins), and the
%                            laminar boundaries of each area.
%   spectrolaminar_data.mat  Spectrolaminar (LFP) data.
%
% Sections (manuscript figures)
%    1. Data preparation
%    2. Prototypical cortex                        Figure 3
%    3. Cross-area conservation                    Figure 2B, Figure S3
%    4. PCA and bootstrap ratios                   Figure 4A-C, Figure S4
%    5. Procrustes analysis                        Figure 4D
%    6. PLSR on individual areas                   Figure S7
%    7. PLSR with leave-one-area-out CV            Figure 5
%    8. PLSR on frequency bands                    Figure S5
%    9. Schizophrenia perturbation                 Figure 6
%
% Requirements: MATLAB R2021a or later, Statistics and Machine Learning Toolbox.

clc; clear; close all;
repo_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_dir, 'functions'));

N_PERMS      = 5000;             % permutations for all permutation tests
N_BOOTSTRAPS = 5000;             % bootstrap resamples for PCA stability
n_components = 4;                % PLSR components

%% 1. Data preparation
slm = load(fullfile(repo_dir, 'spectrolaminar_data.mat'));
pub = load(fullfile(repo_dir, 'LANSA_processed_cyto.mat'));

included_areas = pub.included_areas;     % V1 V3 V4 MT MST LIP TEO PMD 7A 8A PFC
cell_types     = pub.cell_types;         % PV CB CR NG
standardEdges  = pub.standard_edges;     % laminar bin positions (proportion of depth from mid-L4)
anatomy_data   = pub.anatomy_data;       % X: 21 bins x 4 cell types x 11 areas (cell density)
ephys_data     = slm.ephys_data;         % Y: 21 bins x 150 frequencies x 11 areas (relative power)

[num_bins, num_cells, num_areas] = size(anatomy_data);
num_freqs = size(ephys_data, 2);
if ~isequal(size(ephys_data), [num_bins, 150, num_areas])
    error('LANSA:ephysSize', 'ephys_data must be %d x 150 x %d (bins x frequencies x areas).', num_bins, num_areas);
end

% Mean laminar boundaries across areas (layer lines on plots; layer assignment
% in section 9). PMD (agranular) has no L4 boundaries
all_layer_bounds = zeros(7, num_areas);
for ii = 1:num_areas
    if ~strcmp(included_areas{ii}, 'PMD')
        all_layer_bounds(:, ii) = pub.laminar_bounds{ii};
    end
end
mean_layer_bounds = mean(all_layer_bounds, 2);
internal_layers   = mean_layer_bounds(2:end-1);
midpoints         = (mean_layer_bounds(1:end-1) + mean_layer_bounds(2:end)) / 2;
layer_labels      = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'};

% Prototypical cortex: mean across areas
proto_anatomy = mean(anatomy_data, 3);   % 21 x 4
proto_ephys   = mean(ephys_data, 3);     % 21 x 150

% Frequency bands
freq_band_labels = {'Delta', 'Theta', 'Alpha', 'Beta', 'Low Gamma', 'High Gamma'};
bands = {1:3, 4:7, 8:12, 13:30, 31:90, 91:150};
freq_boundaries = [0, 3, 7, 12, 30, 90, 150];

% Colors
cell_colors = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556];   % PV CB CR NG
PC1_COLOR = [0.85 0.325 0.098];   % orange
PC2_COLOR = [0 0.447 0.741];      % blue

%% 2. Prototypical cortex (Figure 3)
rel_anatomy = proto_anatomy ./ max(proto_anatomy);
rel_anatomy_sem2 = 2 * std(anatomy_data, 0, 3) ./ sqrt(num_areas) ./ max(proto_anatomy);   % 2 SEM

figure('Color', 'w');
sgtitle('Mean Anatomy and Electrophysiology');
subplot(1, 4, 1);
hold on;
for i = 1:num_cells
    plotShaded(standardEdges, rel_anatomy(:, i), rel_anatomy_sem2(:, i), cell_colors(i, :), ...
        'DisplayName', cell_types{i});
end
set(gca, 'view', [90 90]);
xline(internal_layers, '--', 'HandleVisibility', 'off');
text(midpoints, -0.05 * ones(size(midpoints)), layer_labels, 'FontWeight', 'bold');
xlim(standardEdges([1 end]));
ylim([-0.1 1.4]);
legend('Location', 'southeast');
title('Mean Relative Density (\pm 2 SEM)');
xlabel('Cortical Position from L4');
ylabel('Relative Density');
styleAxes(gca);

subplot(1, 4, 2:3);
imagesc(1:num_freqs, standardEdges, proto_ephys);
clbr = colorbar;
clbr.Label.String = 'Rel Pow';
clim([min(proto_ephys(:)), max(proto_ephys(:))]);
yline(internal_layers, '--', 'HandleVisibility', 'off');
text(2 * ones(size(midpoints)), midpoints, layer_labels, 'FontWeight', 'bold');
ylim(standardEdges([1 end]));
title('Mean Relative Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Cortical Position from L4');
styleAxes(gca);

gamma      = mean(proto_ephys(:, 31:150), 2);
alpha_beta = mean(proto_ephys(:, 8:30), 2);
subplot(1, 4, 4);
hold on;
plot(standardEdges, gamma, 'r', 'LineWidth', 2, 'DisplayName', '\gamma (31-150 Hz)');
plot(standardEdges, alpha_beta, 'b', 'LineWidth', 2, 'DisplayName', '\alpha\beta (8-30 Hz)');
set(gca, 'view', [90 90]);
xlim(standardEdges([1 end]));
ylim([0.5 1]);
legend('Location', 'best');
title('Mean Relative Power');
ylabel('Relative Power');
styleAxes(gca);

% E/I balance: pyramidal (NG) / interneuron (PV + CB + CR) density
E_I = squeeze(anatomy_data(:, 4, :) ./ sum(anatomy_data(:, 1:3, :), 2));   % 21 x areas
mean_E_I = mean(E_I, 2);
SEM2_E_I = 2 * std(E_I, 0, 2) ./ sqrt(num_areas);

figure('Color', 'w');
plotShaded(standardEdges, mean_E_I / max(mean_E_I), SEM2_E_I / max(mean_E_I), [0.5930 0.4887 0.3213], ...
    'DisplayName', 'E/I');
set(gca, 'view', [90 90]);
xline(internal_layers, '--', 'HandleVisibility', 'off');
text(midpoints, -0.05 * ones(size(midpoints)), layer_labels, 'FontWeight', 'bold');
ylim([-0.1 1.4]);
title('Mean E/I Balance (\pm 2 SEM)');
xlabel('Cortical Position from L4');
ylabel('Relative E/I');
styleAxes(gca);

%% 3. Cross-area conservation (Figure 2B, Figure S3)
% Pearson r between every pair of areas (55 pairs): for each cell type, between
% the 21-bin density profiles; for the spectrolaminar motif, between each area's
% 21 x 150 map unrolled into one vector (also within each frequency band).
% Conservation = mean pairwise r. Significance: permutation test in which each
% area's profile is circularly shifted along depth (permTestMeanR).
rng(42);
summary_names = {};
summary_vals  = [];   % [mean r, median r, p]

figure('Color', 'w');
sgtitle('Conservation of Cell Types Across Cortical Areas');
for jj = 1:num_cells
    [p_val, mean_r, R, median_r] = permTestMeanR(reshape(anatomy_data(:, jj, :), num_bins, 1, num_areas), N_PERMS);
    R(~tril(true(num_areas))) = NaN;
    subplot(2, 2, jj);
    hmap = heatmap(R);
    hmap.Colormap = turbo;
    hmap.ColorLimits = [-1 1];
    hmap.XDisplayLabels = included_areas;
    hmap.YDisplayLabels = included_areas;
    hmap.Title = sprintf('%s | Mean r = %.2f (p = %.3g)', cell_types{jj}, mean_r, p_val);
    summary_names{end+1} = cell_types{jj};
    summary_vals(end+1, :) = [mean_r, median_r, p_val];
end

[p_val, mean_r, R, median_r] = permTestMeanR(ephys_data, N_PERMS);
R(~tril(true(num_areas))) = NaN;
figure('Color', 'w');
hmap = heatmap(R);
hmap.Colormap = turbo;
hmap.ColorLimits = [-1 1];
hmap.XDisplayLabels = included_areas;
hmap.YDisplayLabels = included_areas;
hmap.Title = sprintf('Conservation of SLM | Mean r = %.2f (p = %.3g)', mean_r, p_val);
summary_names{end+1} = 'SLM';
summary_vals(end+1, :) = [mean_r, median_r, p_val];

figure('Color', 'w');
sgtitle('Conservation of SLM Across Cortical Areas (Bands)');
for bb = 1:numel(bands)
    [p_val, mean_r, R, median_r] = permTestMeanR(ephys_data(:, bands{bb}, :), N_PERMS);
    R(~tril(true(num_areas))) = NaN;
    subplot(2, 3, bb);
    hmap = heatmap(R);
    hmap.Colormap = turbo;
    hmap.ColorLimits = [-1 1];
    hmap.XDisplayLabels = included_areas;
    hmap.YDisplayLabels = included_areas;
    hmap.Title = {freq_band_labels{bb}, sprintf('Mean r = %.2f (p = %.3g)', mean_r, p_val)};
    summary_names{end+1} = freq_band_labels{bb}; 
    summary_vals(end+1, :) = [mean_r, median_r, p_val];
end

fprintf('\n--- Conservation across areas ---\n%-12s %8s %9s %10s\n', '', 'Mean r', 'Median r', 'p');
for ii = 1:numel(summary_names)
    fprintf('%-12s %8.2f %9.2f %10.4f\n', summary_names{ii}, summary_vals(ii, :));
end

%% 4. Principal Component Analysis (Figure 4A-C, Figure S4)
% Coefficients (loadings): contribution of each cell type / frequency to a PC.
% Scores: laminar profile of each PC. Explained: % variance per PC.
X_pca = zscore(proto_anatomy);   % z-score each cell type across depth
Y_pca = zscore(proto_ephys);     % z-score each frequency across depth
N_OBS_pca = size(X_pca, 1);
rng(42);

[coeff_X, score_X, ~, ~, explained_X] = pca(X_pca);
[coeff_Y, score_Y, ~, ~, explained_Y] = pca(Y_pca);
fprintf('\n--- PCA on Cell Types (X) Explained Variance ---\n PC1: %.2f%%\n PC2: %.2f%%\n PC3: %.2f%%\n PC4: %.2f%%\n', explained_X(1:4));
fprintf('--- PCA on Frequency Power (Y) Explained Variance ---\n PC1: %.2f%%\n PC2: %.2f%%\n PC3: %.2f%%\n PC4: %.2f%%\n', explained_Y(1:4));

% Retain PCs until cumulative explained variance >= 90% (both datasets: 2 PCs)
explained_variance_threshold = 90;
k_perm = max(find(cumsum(explained_X) >= explained_variance_threshold, 1), ...
             find(cumsum(explained_Y) >= explained_variance_threshold, 1));
fprintf('PCs retained (>= %d%% variance): %d\n', explained_variance_threshold, k_perm);

% --- Figure 4A, 4B: loadings and scores of PC1 and PC2 ---
pc_colors = [PC1_COLOR; PC2_COLOR];
figure('Color', 'w');
sgtitle('Latent Gradients Across Cortical Depth (PCA)');
for pc = 1:2
    subplot(2, 4, pc);
    b = barh(coeff_X(:, pc));
    set(b, 'FaceColor', pc_colors(pc, :), 'EdgeColor', pc_colors(pc, :));
    set(gca, 'YDir', 'reverse', 'YTick', 1:num_cells, 'YTickLabel', cell_types);
    xlim([-1 1]);
    xline(0);
    title(sprintf('PC%d_X (%.1f%%)', pc, explained_X(pc)));
    xlabel('Loading');
    ylabel('Cell Type');
    styleAxes(gca);

    subplot(2, 4, 4 + pc);
    b = barh(coeff_Y(:, pc));
    set(b, 'FaceColor', pc_colors(pc, :), 'EdgeColor', pc_colors(pc, :));
    set(gca, 'YDir', 'reverse', 'YTick', 0:10:num_freqs);
    yline(freq_boundaries, '--');
    xlim([-0.2 0.2]);
    xline(0);
    title(sprintf('PC%d_Y (%.1f%%)', pc, explained_Y(pc)));
    xlabel('Loading');
    ylabel('Frequency (Hz)');
    styleAxes(gca);
end
score_sets = {score_X, score_Y};
score_titles = {'Cell Type Gradient', 'Functional PC Gradients'};
for row = 1:2
    subplot(2, 4, 4 * (row - 1) + (3:4));
    s = score_sets{row};
    p = plot(s(:, 1:2), standardEdges, 'LineWidth', 2);
    set(p(1), 'Color', PC1_COLOR);
    set(p(2), 'Color', PC2_COLOR);
    set(gca, 'YDir', 'reverse');
    ylim(standardEdges([1 end]));
    yline(internal_layers, '--', 'HandleVisibility', 'off');
    text(min(min(s(:, 1:2))) * ones(size(midpoints)), midpoints, layer_labels, 'FontWeight', 'bold');
    xline(0);
    legend(p, {'PC1', 'PC2'}, 'Location', 'southeast');
    title(score_titles{row});
    xlabel('Scores');
    ylabel('Cortical Depth');
    styleAxes(gca);
end

% --- Figure 4C: scree plot ---
figure('Color', 'w');
hold on;
plot(1:numel(explained_X), cumsum(explained_X), '-o', 'Color', [0 0.447 0.741], 'LineWidth', 2);
plot(1:numel(explained_Y), cumsum(explained_Y), '-o', 'Color', [0.8 0.1 0.1], 'LineWidth', 2);
yline(explained_variance_threshold, 'k-');
legend({'Structural Variance', 'Functional Variance'}, 'Location', 'southeast');
title('Scree Plot for X and Y PCAs');
xlabel('Dimensions');
ylabel('Percentage of Explained Variance');
ylim([50 100]);
styleAxes(gca);

% --- Figure S4: bootstrap ratios (BSR) of the loadings ---
% Laminar bins are resampled with replacement; each bootstrap PC is sign-aligned
% to the observed PC. BSR = mean / SD across bootstraps; |BSR| >= 2.5 is stable.
BSR_THRESHOLD = 2.5;
boot_coeff_X = zeros(N_BOOTSTRAPS, size(coeff_X, 1), k_perm);
boot_coeff_Y = zeros(N_BOOTSTRAPS, size(coeff_Y, 1), k_perm);
for i = 1:N_BOOTSTRAPS
    sample_indices = datasample(1:N_OBS_pca, N_OBS_pca, 'Replace', true);
    coeff_boot_X = pca(X_pca(sample_indices, :));
    coeff_boot_Y = pca(Y_pca(sample_indices, :));
    for j = 1:k_perm
        if dot(coeff_X(:, j), coeff_boot_X(:, j)) < 0, coeff_boot_X(:, j) = -coeff_boot_X(:, j); end
        if dot(coeff_Y(:, j), coeff_boot_Y(:, j)) < 0, coeff_boot_Y(:, j) = -coeff_boot_Y(:, j); end
        boot_coeff_X(i, :, j) = coeff_boot_X(:, j);
        boot_coeff_Y(i, :, j) = coeff_boot_Y(:, j);
    end
end
BSR_X = squeeze(mean(boot_coeff_X, 1) ./ std(boot_coeff_X, 0, 1));   % cell types x PCs
BSR_Y = squeeze(mean(boot_coeff_Y, 1) ./ std(boot_coeff_Y, 0, 1));   % frequencies x PCs

figure('Color', 'w');
sgtitle('PCA Bootstrap Ratio for Coefficients');
for pc = 1:2
    subplot(2, 2, pc);
    b = barh(BSR_X(:, pc));
    set(b, 'FaceColor', pc_colors(pc, :), 'EdgeColor', pc_colors(pc, :));
    set(gca, 'YDir', 'reverse', 'YTick', 1:num_cells, 'YTickLabel', cell_types);
    xline([-BSR_THRESHOLD, BSR_THRESHOLD], 'r--', 'LineWidth', 1.5);
    xline(0);
    title(sprintf('PC%d_X BSR', pc));
    xlabel('Bootstrap Ratio');
    ylabel('Cell Type');
    styleAxes(gca);

    subplot(2, 2, 2 + pc);
    b = barh(BSR_Y(:, pc));
    set(b, 'FaceColor', pc_colors(pc, :), 'EdgeColor', pc_colors(pc, :));
    set(gca, 'YDir', 'reverse', 'YTick', 0:10:num_freqs);
    yline(freq_boundaries, '--');
    xline([-BSR_THRESHOLD, BSR_THRESHOLD], 'r--', 'LineWidth', 1.5);
    xline(0);
    title(sprintf('PC%d_Y BSR', pc));
    xlabel('Bootstrap Ratio');
    ylabel('Frequency (Hz)');
    styleAxes(gca);
end

%% 5. Procrustes Analysis (Figure 4D)
% PCA on each area separately; Procrustes distance d between the area's PC1-2
% scores and the global PC1-2 scores (0 = identical, 1 = dissimilar).
% Null: laminar bins of each area permuted N_PERMS times; critical d = 1st percentile.
d_scores = zeros(num_areas, 2);   % [structure (X), function (Y)]
d_perm_X = zeros(num_areas, N_PERMS);
d_perm_Y = zeros(num_areas, N_PERMS);
for t = 1:num_areas
    X_single_area = zscore(anatomy_data(:, :, t));
    Y_single_area = zscore(ephys_data(:, :, t));
    [~, score_anatomy] = pca(X_single_area);
    [~, score_ephys]   = pca(Y_single_area);
    d_scores(t, 1) = procrustes(score_X(:, 1:k_perm), score_anatomy(:, 1:k_perm));
    d_scores(t, 2) = procrustes(score_Y(:, 1:k_perm), score_ephys(:, 1:k_perm));
    for p = 1:N_PERMS
        p_indices = randperm(N_OBS_pca);
        [~, score_perm_X] = pca(X_single_area(p_indices, :));
        [~, score_perm_Y] = pca(Y_single_area(p_indices, :));
        d_perm_X(t, p) = procrustes(score_X(:, 1:k_perm), score_perm_X(:, 1:k_perm));
        d_perm_Y(t, p) = procrustes(score_Y(:, 1:k_perm), score_perm_Y(:, 1:k_perm));
    end
end
d_critical = [prctile(d_perm_X(:), 1), prctile(d_perm_Y(:), 1)];

fprintf('\n--- Procrustes distance to global PC1-2 ---\n');
for t = 1:num_areas
    fprintf('%-4s d(X) = %.3f  d(Y) = %.3f\n', included_areas{t}, d_scores(t, 1), d_scores(t, 2));
end
fprintf('Critical d (1st percentile of null): d(X) = %.3f, d(Y) = %.3f\n', d_critical);

figure('Color', 'w');
hold on;
scatter(d_scores(:, 1), d_scores(:, 2), 40, 'k', 'filled');
text(d_scores(:, 1) + 0.01, d_scores(:, 2), included_areas, 'FontWeight', 'bold');
quadrant_color = [0.5 0.5 0.5];
text(0.45, 0.45, 'X and Y congruent', 'Color', quadrant_color, 'FontWeight', 'bold');
text(0.45, 0.85, 'X congruent, Y incongruent', 'Color', quadrant_color, 'FontWeight', 'bold');
text(0.85, 0.45, {'X incongruent,', 'Y congruent'}, 'Color', quadrant_color, 'FontWeight', 'bold');
text(0.85, 0.85, {'X and Y', 'incongruent'}, 'Color', quadrant_color, 'FontWeight', 'bold');
xline(d_critical(1), 'r--');
yline(d_critical(2), 'r--');
xlim([0 1]);
ylim([0 1]);
title('Procrustes Distance from Global PC1 and PC2');
xlabel('d (Structure)');
ylabel('d (Function)');
styleAxes(gca);

%% 6. PLSR on Individual Areas (Figure S7)
% Predictors (buildFeatureMatrix): z-scored PV, CB, CR, NG, their 6 pairwise
% interactions and E/I balance (11 predictors). One model per area (response:
% that area's relative power). Each area's model predicts from its own
% cytoarchitecture; the prediction is correlated (2-D Pearson r over all bins
% and frequencies) with every area's observed motif. The beta matrices are
% also correlated across areas.
predictor_labels = {'PV', 'CB', 'CR', 'NG', 'PV×CB', 'PV×CR', 'PV×NG', 'CB×CR', 'CB×NG', 'CR×NG', 'E/I'};
P = numel(predictor_labels);
Q = num_freqs;
Y_FREQ_Hz = 1:Q;
C1 = cell_colors(1, :); C2 = cell_colors(2, :); C3 = cell_colors(3, :); C4 = cell_colors(4, :);
predictor_colors = [C1; C2; C3; C4; (C1+C2)/2; (C1+C3)/2; (C1+C4)/2; (C2+C3)/2; (C2+C4)/2; (C3+C4)/2; (C1+C2+C3)/3];

X_all = cell(num_areas, 1);   % predictor matrix of each area
beta_by_area = zeros(P + 1, Q, num_areas);
r_pred_cross_area = zeros(num_areas);   % rows: predicting area; columns: observed area
for a = 1:num_areas
    X_all{a} = buildFeatureMatrix(anatomy_data(:, :, a));
    [~, ~, ~, ~, BETA_cross] = plsregress(X_all{a}, ephys_data(:, :, a), n_components);
    beta_by_area(:, :, a) = BETA_cross;
    ephys_indv_pred = X_all{a} * BETA_cross(2:end, :);
    for jj = 1:num_areas
        r_pred_cross_area(a, jj) = corr(reshape(ephys_data(:, :, jj), [], 1), ephys_indv_pred(:));
    end
end
r_beta_across_areas = zeros(num_areas);
for a = 1:num_areas
    for jj = 1:num_areas
        r_beta_across_areas(a, jj) = corr(reshape(beta_by_area(:, :, a), [], 1), reshape(beta_by_area(:, :, jj), [], 1));
    end
end
offMask = ~eye(num_areas);
lowMask = tril(true(num_areas), -1);
on_cross   = mean(diag(r_pred_cross_area));
off_cross  = mean(r_pred_cross_area(offMask));
r_beta_val = mean(r_beta_across_areas(lowMask));

figure('Color', 'w');
hmap = heatmap(r_pred_cross_area);
hmap.Colormap = turbo;
hmap.ColorLimits = [-1 1];
hmap.XDisplayLabels = included_areas;
hmap.YDisplayLabels = included_areas;
hmap.Title = {'Cross Area Predictability', sprintf('Mean Internal Predictive Strength r = %.2f', on_cross), ...
    sprintf('Mean External Predictive Strength r = %.2f', off_cross)};

R = r_beta_across_areas;
R(~tril(true(num_areas))) = NaN;
figure('Color', 'w');
hmap = heatmap(R);
hmap.Colormap = turbo;
hmap.ColorLimits = [0.75 1];
hmap.XDisplayLabels = included_areas;
hmap.YDisplayLabels = included_areas;
hmap.Title = sprintf('Conservation of \\beta Weights (mean r = %.3f)', r_beta_val);

% Area-specific beta coefficients
figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
ax = gobjects(1, num_areas);
for jj = 1:num_areas
    ax(jj) = subplot(3, 4, jj);
    hold on;
    yline(0, 'k', 'HandleVisibility', 'off');
    xline(freq_boundaries, 'k:', 'HandleVisibility', 'off');
    for i = 1:P
        line_style = '-';
        if i > 4, line_style = '--'; end
        plot(Y_FREQ_Hz, beta_by_area(i + 1, :, jj), line_style, 'LineWidth', 2, 'Color', predictor_colors(i, :));
    end
    title(included_areas{jj});
    xlim([0 Q]);
    if mod(jj - 1, 4) == 0, ylabel('\beta Coefficient'); end
    if jj > 8, xlabel('Freq (Hz)'); end
    styleAxes(gca);
end
legend_ax = subplot(3, 4, 12);
axis(legend_ax, 'off');
lgd = legend(ax(1), predictor_labels, 'Location', 'none');
set(lgd, 'Units', 'normalized', 'Position', get(legend_ax, 'Position'), 'Box', 'off');
sgtitle('Area-Specific Beta Coefficients');

% --- Permutation test: internal prediction, external prediction, beta conservation ---
% The same random shuffle of laminar bins is applied to every area's predictors.
null_diag = zeros(N_PERMS, num_areas);
null_on   = zeros(N_PERMS, 1);
null_off  = zeros(N_PERMS, 1);
null_beta = zeros(N_PERMS, 1);
for p = 1:N_PERMS
    idx = randperm(num_bins);
    r_perm = zeros(num_areas);
    B_perm = zeros(size(beta_by_area));
    for a = 1:num_areas
        Xp = X_all{a}(idx, :);
        [~, ~, ~, ~, B] = plsregress(Xp, ephys_data(:, :, a), n_components);
        B_perm(:, :, a) = B;
        pred = Xp * B(2:end, :);
        for jj = 1:num_areas
            r_perm(a, jj) = corr(reshape(ephys_data(:, :, jj), [], 1), pred(:));
        end
    end
    rb = zeros(num_areas);
    for i = 1:num_areas
        for j = 1:i-1
            rb(i, j) = corr(reshape(B_perm(:, :, i), [], 1), reshape(B_perm(:, :, j), [], 1));
        end
    end
    null_diag(p, :) = diag(r_perm)';
    null_on(p)      = mean(diag(r_perm));
    null_off(p)     = mean(r_perm(offMask));
    null_beta(p)    = mean(rb(lowMask));
end
permP = @(null_dist, obs) (sum(null_dist >= obs) + 1) / (N_PERMS + 1);

fprintf('\n--- PLSR on individual areas ---\n');
for a = 1:num_areas
    fprintf('Area %s: Observed r = %.4f, p = %.4f\n', included_areas{a}, r_pred_cross_area(a, a), ...
        permP(null_diag(:, a), r_pred_cross_area(a, a)));
end
fprintf('Internal prediction: mean r = %.4f, p = %.4f\n', on_cross, permP(null_on, on_cross));
fprintf('External prediction: mean r = %.4f, p = %.4f\n', off_cross, permP(null_off, off_cross));
fprintf('Beta conservation:   mean r = %.4f, p = %.4f\n', r_beta_val, permP(null_beta, r_beta_val));

%% 7. PLSR with Leave-One-Area-Out Cross Validation (Figure 5)
% Training: mean cytoarchitecture and mean (z-scored) motif of the other 10
% areas. Test: the held-out area. Q^2 = 1 - PRESS / TSS, with TSS around the
% training mean. Q^2 significance per area and frequency: training motif rows
% shuffled N_PERMS times; p-values combined with Fisher's method.
X_plsr = buildFeatureMatrix(proto_anatomy);

beta_matrix_loo    = zeros(P + 1, Q, num_areas);
R_squared_internal = zeros(num_areas, Q);   % Q^2
r_prediction       = zeros(num_areas, Q);
p_val_plsr         = zeros(num_areas, Q);
XL_loo = zeros(P, n_components, num_areas);   XS_loo = zeros(num_bins, n_components, num_areas);
YL_loo = zeros(Q, n_components, num_areas);   YS_loo = zeros(num_bins, n_components, num_areas);
PCTVAR_loo = zeros(2, n_components, num_areas);

for t = 1:num_areas
    train_idx = setdiff(1:num_areas, t);
    X_train_final = buildFeatureMatrix(mean(anatomy_data(:, :, train_idx), 3));
    Y_train_z     = zscore(mean(ephys_data(:, :, train_idx), 3));
    X_test_final  = buildFeatureMatrix(anatomy_data(:, :, t));
    Y_test_z      = zscore(ephys_data(:, :, t));

    [XL, YL, XS, YS, BETA, PCTVAR] = plsregress(X_train_final, Y_train_z, n_components);
    beta_matrix_loo(:, :, t) = BETA;
    XL_loo(:, :, t) = XL;  YL_loo(:, :, t) = YL;
    XS_loo(:, :, t) = XS;  YS_loo(:, :, t) = YS;
    PCTVAR_loo(:, :, t) = PCTVAR;

    Y_test_pred = BETA(1, :) + X_test_final * BETA(2:end, :);
    PRESS = sum((Y_test_z - Y_test_pred).^2, 1);
    TSS   = sum((Y_test_z - mean(Y_train_z, 1)).^2, 1);
    R_squared_internal(t, :) = 1 - PRESS ./ TSS;

    R_perm_set = zeros(N_PERMS, Q);
    for n = 1:N_PERMS
        Y_shuff = Y_train_z(randperm(num_bins), :);
        [~, ~, ~, ~, B_perm] = plsregress(X_train_final, Y_shuff, n_components);
        Y_shuff_pred = B_perm(1, :) + X_test_final * B_perm(2:end, :);
        R_perm_set(n, :) = 1 - sum((Y_test_z - Y_shuff_pred).^2, 1) ./ TSS;
    end

    for f = 1:Q
        p_val_plsr(t, f) = sum(R_perm_set(:, f) >= R_squared_internal(t, f)) / N_PERMS;
        r_prediction(t, f) = corr(Y_test_z(:, f), Y_test_pred(:, f));
    end
end

% Fisher's combined probability test
all_p = max(p_val_plsr(:), realmin);
overall_p_Fisher = 1 - chi2cdf(-2 * sum(log(all_p)), 2 * numel(all_p));

mean_pearson_r = mean(r_prediction, 1);
mean_q2 = mean(R_squared_internal, 1);
mean_PCTVAR = mean(PCTVAR_loo, 3) * 100;
fprintf('\n--- PLSR LOOCV (%d components) ---\n', n_components);
fprintf('Mean r = %.2f, mean Q^2 = %.2f, Fisher combined p = %.3g\n', mean(mean_pearson_r), mean(mean_q2), overall_p_Fisher);
fprintf('Variance explained: X = %.2f%%, Y = %.2f%%\n', sum(mean_PCTVAR(1, :)), sum(mean_PCTVAR(2, :)));

% Band summary: r and Q^2 averaged within band per held-out area; two-sided
% Wilcoxon signed-rank test against zero across areas
band_labels_loo = {'Delta', 'Theta', 'Alpha_Beta', 'Lo_Gamma', 'Hi_Gamma'};
freq_bands_loo  = {1:3, 4:7, 8:30, 31:90, 91:150};
fprintf('%-11s %14s %9s %14s %9s\n', 'Band', 'Q^2 range', 'p(Q^2)', 'r range', 'p(r)');
for b = 1:numel(freq_bands_loo)
    fb = freq_bands_loo{b};
    fprintf('%-11s %6.2f to %5.2f %9.4f %6.2f to %5.2f %9.4f\n', band_labels_loo{b}, ...
        min(mean_q2(fb)), max(mean_q2(fb)), signrank(mean(R_squared_internal(:, fb), 2)), ...
        min(mean_pearson_r(fb)), max(mean_pearson_r(fb)), signrank(mean(r_prediction(:, fb), 2)));
end

% --- Figure 5A: beta coefficients (mean +/- 2 SE across LOOCV folds) ---
beta_mean = mean(beta_matrix_loo, 3);
beta_variance = 2 * std(beta_matrix_loo, 0, 3) ./ sqrt(num_areas);
figure('Color', 'w');
hold on;
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(freq_boundaries, 'k--', 'HandleVisibility', 'off');
for i = 1:P
    h = plotShaded(Y_FREQ_Hz, beta_mean(i + 1, :), beta_variance(i + 1, :), predictor_colors(i, :), ...
        'DisplayName', predictor_labels{i});
    if i > 4, set(h, 'LineStyle', '--'); end
end
set(findobj(gca, 'Type', 'patch'), 'FaceAlpha', 0.15);
ylim([-0.5 0.5]);
legend('Location', 'eastoutside');
title({'Predictive Strength of Cell Types and Interactions Across Frequencies', sprintf('Components = %d', n_components)});
xlabel('Frequency (Hz)');
ylabel('\beta Coefficient (Mean \pm 2SE)');
styleAxes(gca);

% --- Figure 5B (Pearson r) and 5C (Q^2) across frequencies ---
metric_data   = {r_prediction, R_squared_internal};
metric_ylabel = {'Pearson Correlation (r \pm 2SE)', 'Internal Predictive Strength (Q^2 \pm 2SE)'};
metric_title  = {'Pearson Correlation Between Observed and Predicted Y', 'Internal Predictive Strength of Model'};
metric_name   = {'r', 'Q^2'};
for m = 1:2
    vals = metric_data{m};
    mean_curve = mean(vals, 1);
    figure('Color', 'w');
    hold on;
    plotShaded(Y_FREQ_Hz, mean_curve, 2 * std(vals, 0, 1) ./ sqrt(num_areas), [0 0 0], 'HandleVisibility', 'off');
    for a = 1:num_areas
        line_style = '-';
        if a > 6, line_style = '--'; end
        plot(Y_FREQ_Hz, vals(a, :), line_style, 'LineWidth', 1, 'DisplayName', included_areas{a});
    end
    xline(freq_boundaries, '-', 'Color', [0.83 0.83 0.83], 'HandleVisibility', 'off');
    yline(0, '-', 'Color', [0.83 0.83 0.83], 'HandleVisibility', 'off');
    yline(mean(mean_curve), 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
    text(Q, mean(mean_curve), sprintf('Mean %s = %.2f', metric_name{m}, mean(mean_curve)), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
    if m == 1, ylim([-1.25 1.25]); end
    legend('Location', 'eastoutside');
    title({metric_title{m}, sprintf('Leave-One-Area-Out Cross Validation (N COMPS = %d)', n_components)});
    xlabel('Frequency (Hz)');
    ylabel(metric_ylabel{m});
    styleAxes(gca);
end

% --- Predicted spectrolaminar motif (prototypical cortex, mean LOOCV beta) ---
Y_pred_avg_raw = (X_plsr * beta_mean(2:end, :) + beta_mean(1, :)) .* std(ephys_data, 0, 3) + mean(ephys_data, 3);
figure('Color', 'w');
imagesc(Y_FREQ_Hz, standardEdges, Y_pred_avg_raw);
clbr = colorbar;
clbr.Label.String = 'Rel Pow';
clim([0 1]);
title({'Predicted Electrophysiology', 'Reversing z-score Transformation'});
xlabel('Frequency (Hz)');
ylabel('Cortical Depth (Proportion from L4)');
styleAxes(gca);

% --- PLSR model internal architecture (mean +/- 2 SE across folds) ---
comp_labels = cell(1, n_components);
for jj = 1:n_components
    comp_labels{jj} = sprintf('Comp %d: PEV_X = %0.1f%%, PEV_Y = %0.1f%%', jj, mean_PCTVAR(1, jj), mean_PCTVAR(2, jj));
end
arch_data   = {XL_loo, XS_loo, YL_loo, YS_loo};
arch_x      = {1:P, standardEdges, Y_FREQ_Hz, standardEdges};
arch_titles = {'Predictor Weights (XL)', 'Observation Scores (XS)', 'Response Weights (YL)', 'Observation Scores (YS)'};
figure('Color', 'w');
sgtitle('PLSR Model Internal Architecture');
for k = 1:4
    subplot(2, 2, k);
    hold on;
    mu = mean(arch_data{k}, 3);
    se2 = 2 * std(arch_data{k}, 0, 3) ./ sqrt(num_areas);
    for jj = 1:n_components
        plotShaded(arch_x{k}, mu(:, jj), se2(:, jj), predictor_colors(jj, :), 'LineWidth', 1.5);
    end
    if k == 1
        set(gca, 'XTick', 1:P, 'XTickLabel', predictor_labels);
        xtickangle(45);
        ylabel('X-Loadings');
    elseif k == 3
        xlabel('Frequency (Hz)');
        ylabel('Y-Loadings');
        legend(comp_labels, 'Location', 'best');
    else
        set(gca, 'view', [90 90]);
        xline(internal_layers, '--', 'HandleVisibility', 'off');
        xlabel('Cortical Depth');
        ylabel('Contribution');
    end
    yline(0, 'k--', 'HandleVisibility', 'off');
    title(arch_titles{k});
    styleAxes(gca);
end

%% 8. PLSR Frequency Band Analysis (Figure S5)
% LOOCV repeated with z-scored power averaged within each of six bands.
nBands = numel(bands);
beta_matrix_bands = zeros(num_areas, P + 1, nBands);
for t = 1:num_areas
    loo_idx = setdiff(1:num_areas, t);
    X_band = buildFeatureMatrix(mean(anatomy_data(:, :, loo_idx), 3));
    Y_freqs_z = zscore(mean(ephys_data(:, :, loo_idx), 3));
    Y_band = zeros(num_bins, nBands);
    for b = 1:nBands
        Y_band(:, b) = mean(Y_freqs_z(:, bands{b}), 2);
    end
    [~, ~, ~, ~, BETA_band] = plsregress(X_band, Y_band, n_components);
    beta_matrix_bands(t, :, :) = BETA_band;
end
beta_mean_bands = squeeze(mean(beta_matrix_bands(:, 2:end, :), 1));                      % predictors x bands
beta_variance_bands = 2 * squeeze(std(beta_matrix_bands(:, 2:end, :), 0, 1)) ./ sqrt(num_areas);

figure('Color', 'w');
hold on;
b = bar(beta_mean_bands');
for bb = 1:P
    set(b(bb), 'FaceColor', predictor_colors(bb, :), 'BarWidth', 0.95);
    errorbar(b(bb).XEndPoints, beta_mean_bands(bb, :), beta_variance_bands(bb, :), 'k', ...
        'LineStyle', 'none', 'CapSize', 4, 'LineWidth', 2, 'HandleVisibility', 'off');
end
yline([-0.1 0.1], 'r', 'HandleVisibility', 'off');
legend(b, predictor_labels, 'Location', 'eastoutside');
set(gca, 'XTick', 1:nBands, 'XTickLabel', freq_band_labels);
ylim([-0.5 0.5]);
title('Mean Beta Weight Across Freq Bands');
xlabel('Frequency Band');
ylabel('Beta Weight \pm 2 SEM (Variance Across Areas)');
styleAxes(gca);

%% 9. Schizophrenia Perturbation (Figure 6)
% PV, CB and CR densities are knocked down layer by layer by the Hedges' g
% effect sizes of Mulvey et al. (2026): factor = (mean + g * SD) / mean, with
% the mean and SD of each cell type pooled across all bins and areas. Each
% area's control and knocked-down cytoarchitecture are passed through its
% area-specific model (beta_by_area, section 6); predictions are returned to
% relative-power units with the across-area mean and SD of the motif.
% Significance: two-sided Wilcoxon signed-rank test across areas, per frequency.
pooled_anatomy = reshape(permute(anatomy_data, [1 3 2]), [], num_cells);
ref.mu    = mean(pooled_anatomy, 1);
ref.sigma = std(pooled_anatomy, 0, 1);
ei_all    = (pooled_anatomy(:, 4) + eps) ./ sum(pooled_anatomy(:, 1:3), 2);
ref.mu_ei    = mean(ei_all);
ref.sigma_ei = std(ei_all);

mu_ephys    = mean(ephys_data, 3);
sigma_ephys = std(ephys_data, 0, 3);

% Hedges' g from the meta-analysis, layers L1-L6
g_pv = [-0.11, -0.49, -0.429, -0.24, -0.401, -0.145];
g_cb = [0.161, -0.301, -0.044, 0.027, -0.050, 0.031];
g_cr = [-0.40, -0.028, -0.054, 0.087, 0.078, -0.057];
kd_layer = [(ref.mu(1) + g_pv * ref.sigma(1)) / ref.mu(1);
            (ref.mu(2) + g_cb * ref.sigma(2)) / ref.mu(2);
            (ref.mu(3) + g_cr * ref.sigma(3)) / ref.mu(3)]';   % 6 layers x 3 cell types

% Assign each laminar bin to a layer using the mean laminar boundaries
layer_bounds = [standardEdges(1); internal_layers; standardEdges(end)];
kd_vecs_combined = ones(num_bins, num_cells);   % NG unchanged
for i = 1:6
    layer_mask = (standardEdges >= layer_bounds(i)) & (standardEdges < layer_bounds(i + 1));
    kd_vecs_combined(layer_mask, 1:3) = repmat(kd_layer(i, :), sum(layer_mask), 1);
end
kd_vecs_combined(end, 1:3) = kd_layer(end, :);

perc_change_by_area = zeros(num_areas, Q);
for a = 1:num_areas
    beta_weight = beta_by_area(:, :, a);
    anatomy_area_raw = anatomy_data(:, :, a);
    anatomy_knocked_raw = anatomy_area_raw .* kd_vecs_combined;

    Y_pred_reg_raw = (beta_weight(1, :) + buildFeatureMatrix(anatomy_area_raw, ref) * beta_weight(2:end, :)) .* sigma_ephys + mu_ephys;
    Y_pred_kd_raw  = (beta_weight(1, :) + buildFeatureMatrix(anatomy_knocked_raw, ref) * beta_weight(2:end, :)) .* sigma_ephys + mu_ephys;

    mean_reg = mean(Y_pred_reg_raw, 1);
    perc_change_by_area(a, :) = (mean(Y_pred_kd_raw, 1) - mean_reg) ./ mean_reg * 100;
end
p_freq_across_areas = zeros(1, Q);
for f = 1:Q
    p_freq_across_areas(f) = signrank(perc_change_by_area(:, f));
end
sig_mask_freq_wise = p_freq_across_areas < 0.05;
fprintf('\n--- Schizophrenia perturbation ---\nSignificant frequencies (p < 0.05): %s\n', mat2str(Y_FREQ_Hz(sig_mask_freq_wise)));
[p_min_beta, k_min] = min(p_freq_across_areas(13:30));
fprintf('Beta band (13-30 Hz): smallest p = %.4f at %d Hz\n', p_min_beta, 12 + k_min);

% --- Figure 6A: control vs knocked-down density (prototypical cortex) ---
rel_proto = proto_anatomy ./ max(proto_anatomy);
kd_anatomy = proto_anatomy .* kd_vecs_combined;
rel_kd = kd_anatomy ./ max(proto_anatomy);
figure('Color', 'w');
hold on;
for i = 1:num_cells
    plot(standardEdges, rel_proto(:, i), 'Color', cell_colors(i, :), 'LineWidth', 2, 'DisplayName', [cell_types{i} ' Normal']);
    plot(standardEdges, rel_kd(:, i), '--', 'Color', cell_colors(i, :), 'LineWidth', 2, 'DisplayName', [cell_types{i} ' SZ']);
end
set(gca, 'view', [90 90]);
xline(internal_layers, '--', 'HandleVisibility', 'off');
text(midpoints, 0.05 * ones(size(midpoints)), layer_labels, 'FontWeight', 'bold');
ylim([0 1.1]);
legend('Location', 'southeast');
title('Normal vs. Knocked-Down Cell Density');
xlabel('Cortical Position from L4');
ylabel('Relative Density');
styleAxes(gca);

% --- Figure 6B: frequency-wise percent change in power (mean +/- SEM across areas) ---
SZ_COLOR = [0.38, 0.49, 0.27];
mean_curve = mean(perc_change_by_area, 1);
sem_curve  = std(perc_change_by_area, 0, 1) / sqrt(num_areas);
figure('Color', 'w');
hold on;
plotShaded(Y_FREQ_Hz, mean_curve, sem_curve, SZ_COLOR, 'LineWidth', 2.5, 'DisplayName', 'All Cells');
plot_sig_bars = nan(1, Q);
plot_sig_bars(sig_mask_freq_wise) = 0.29;
plot(Y_FREQ_Hz, plot_sig_bars, 'Color', SZ_COLOR, 'LineWidth', 5, 'HandleVisibility', 'off');
yline(0, '-k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xline(freq_boundaries(1:end-1), '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
xlim([0 Q]);
ylim([-0.3 0.3]);
set(gca, 'XTick', 0:10:Q);
legend('Location', 'southwest', 'Box', 'off');
title('Frequency-Wise Change in Power');
xlabel('Frequency (Hz)');
ylabel('% Change in Power');
styleAxes(gca);
