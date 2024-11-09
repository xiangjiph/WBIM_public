clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'Ablation_test';
exp_name = 'WBIMGel_20240306';
data_root_folder = fullfile(DataManager.fp_layer(exp_group, exp_name, 1), ...
    'Explore');
si_im_cvrt = @(x) 2 * uint16(x);

vis_root_folder = DataManager.fp_visualization_folder(exp_group, exp_name);
vis_file_prefix = sprintf('%s_%s', exp_group, exp_name);
%% Bubble dynamics 117 & 118
% single pulse ablation, 1 J/cm2, frame rate 21.69 Hz. 700 x 700 pixel,
% 0.25 um / pixel). The image is taken about 9.5 seconds after the
% ablation.
% The width of the damage site is ~ 0.5 um at about 30 second after
% ablation. The initial bubble size is about 6 um in width. 
fp_idx = 118;
im_fn = sprintf('WBIM_%05d.tif', fp_idx);
[~, im_name, ~] = fileparts(im_fn);
im_fp = fullfile(data_root_folder, 'Bubble_dynamics', im_fn);

pixel_size = 0.25;
frame_rate_Hz = 21.69;
im = si_im_cvrt(DataManager.load_single_tiff(im_fp));
im_size = size(im);
im_duration = im_size(3) / frame_rate_Hz;

% Crop the relavent region: 
x_range = round(8 / pixel_size * [-1, 1] + im_size(2) / 2);
y_range = round(40 / pixel_size * [-1, 1] + im_size(1) / 2);
% Remove the first frame. Different laser intensity? 
im_c = im(y_range(1) : y_range(2) - 1, x_range(1) : x_range(2) - 1, :);
% The sample vibration is pretty obvious here. 
    %%
vis_frame_idx = round([0.3, 2.3, 4.3, 6.3, 8.3, 10.3, 12.3, 24.3] * frame_rate_Hz);
num_vis = numel(vis_frame_idx);

fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
for i_vis = 1 : num_vis
    ax_hdl = subplot(1, num_vis, i_vis);
    imagesc(im_c(:, :, vis_frame_idx(i_vis)));
    ax_hdl.DataAspectRatio = [1,1,1];
    ax_hdl.Title.String = sprintf('%d s', round(vis_frame_idx(i_vis) / frame_rate_Hz + 9.7));
    
    ax_hdl.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_hdl.XTick, 'UniformOutput', false);
    ax_hdl.YTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_hdl.YTick, 'UniformOutput', false);
    colormap(ax_hdl, 'gray');
    if i_vis == 1
        ax_hdl.XLabel.String = 'X (\mum)';
        ax_hdl.YLabel.String = 'Y (\mum)';
    end
end
vis_fp = fullfile(vis_root_folder, 'Bubble_dynamics', sprintf('%s_%s_bubble_dynamics.png', ...
    vis_file_prefix, im_name));
fun_print_image_in_several_formats(fig_hdl, vis_fp);
%% Ablation volume 
site_idx = 36;
site_folder = sprintf('Site_%d', site_idx);
exp_folder = fullfile(data_root_folder, site_folder);
exp_files = dir(fullfile(exp_folder, '*.tif'));
% Load all files 
exp_im = arrayfun(@(x) si_im_cvrt(DataManager.load_single_tiff(fullfile(exp_folder, ...
    x.name))), exp_files, 'UniformOutput', false);
pixel_size = 1/6;
frame_avg = 5;
pulse_spacing_um = 20;

v_slit_mm = 1.27;
v_slit_fwhm_y_um = v_slit_mm / (700 / 7.2) * 1000;

est_fwhm_z_um = 5;
est_fwhm_y_um = 75;
est_fwhm_x_um = 4;

est_fwhm_y_pxl = round(est_fwhm_y_um / pixel_size);
est_fwhm_x_pxl = round(est_fwhm_x_um / pixel_size);
est_fwhm_z_pxl = round(est_fwhm_z_um);
ctr_y = round(im_size(1) / 2);

registration_Q = true;
%%
% Smooth the image using median filter along the y axis
exp_im_sm = cellfun(@(x) medfilt3(x, [5, 1, 1]), exp_im, 'UniformOutput', false);
exp_im_avg = cellfun(@(x) fun_image_stack_average_frame(x, frame_avg, @mean), ...
    exp_im_sm, 'UniformOutput', false);
% Move the tile after ablation up? 
%% Registration - for potential global deformation 
% No global translation detected.
fix_idx = 1;
reg_idx = 3;
if registration_Q
    tmp_ref = exp_im_avg{fix_idx};
    tmp_reg = exp_im_avg{reg_idx};
    tmp_ref_mask = tmp_ref >= 1e4;
    tmp_reg_mask = tmp_reg >= 1e4;
    fft_search_range = [10, 10, 10];
    [disp_vec, tmp2] = MaskedTranslationRegistration(tmp_ref, tmp_reg, ...
        tmp_ref_mask, tmp_reg_mask, fft_search_range);
    if any(disp_vec)
        tmp_reg_mv = imtranslate(tmp_reg, disp_vec);
        tmp_empty_mask = (tmp_reg_mv == 0);
        tmp_reg_mv(tmp_empty_mask) = tmp_reg(tmp_empty_mask);
        
        % Check the zy profile
        fig_hdl = figure;
        ax_1 = subplot(1,2,1);
        imshowpair(squeeze(tmp_ref(:, 300, :)), squeeze(tmp_reg(:, 300, :)));
        ax_2 = subplot(1,2,2);
        imshowpair(squeeze(tmp_ref(:, 300, :)), squeeze(tmp_reg_mv(:, 300, :)));
        exp_im_avg{reg_idx} = tmp_reg_mv;
    end
end
im_diff = single(exp_im_avg{1}) - single(exp_im_avg{reg_idx});
im_size = size(im_diff);
%% Compute z-distance between the ablation site to the surface of the gelatin
% The image (site 40) has two peaks in the histogram - well-suited for otsu
bg_mask = ~imbinarize(exp_im_avg{1});
max_bg_depth = size(bg_mask, 3) - 1;
bg_mask = imclose(bg_mask, strel('sphere', 3));
% bg_mask_cc = bwconncomp(bg_mask);
% [bg_mask_cc_size, cc_ind] = sort(cellfun(@numel, bg_mask_cc.PixelIdxList), 'descend');
% bg_mask(cat(1, bg_mask_cc.PixelIdxList{cc_ind(bg_mask_cc_size < 2e4)})) = false;
% bg_mask = imclose(bg_mask, strel('cuboid', [25, 25, 2]));
bg_mask_zcount = sum(bg_mask, 3);

% Find the boundary voxels in z direction
bg_mask_g = (bg_mask(:, :, 2:end) - bg_mask(:, :, 1:end-1)) < 0;
% bg_mask_sm_g = padarray(bg_mask_sm_g, [0, 0, 1], 'post');
bg_mask_edge_ind = find(bg_mask_g(:, :, 1:max_bg_depth));
bg_mask_edge_sub = fun_ind2sub(size(bg_mask_g), bg_mask_edge_ind);

bg_mask_edge_sub_m = mean(bg_mask_edge_sub, 1);
% Assume z is a function of x and y
X_mat = cat(2, bg_mask_edge_sub(:, [1,2]), ones(numel(bg_mask_edge_ind), 1));
n_vec = inv(X_mat.' * X_mat) * (X_mat.' * bg_mask_edge_sub(:, 3));
bg_mask_sm_p = bg_mask;
bg_mask_sm_p_ind = find(bg_mask_sm_p);
bg_mask_sm_p_sub = fun_ind2sub(im_size, bg_mask_sm_p_ind);
is_outlier_Q = (bg_mask_sm_p_sub(:, 3) - bg_mask_sm_p_sub(:, 1) .* n_vec(1) - ...
    bg_mask_sm_p_sub(:, 2) .* n_vec(2) - n_vec(3) ) > 2;
bg_mask_sm_p(bg_mask_sm_p_ind(is_outlier_Q)) = false;

[plt_x, plt_y] = meshgrid(1:im_size(2), 1:im_size(1));
plt_z = plt_x .* n_vec(2) + plt_y .* n_vec(1) + n_vec(3);
% Visualization
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
surf(ax_hdl, plt_x, plt_y, plt_z);
hold(ax_hdl, 'on');
scatter3(ax_hdl, bg_mask_edge_sub(1:100:end, 2), bg_mask_edge_sub(1:100:end, 1), ...
    bg_mask_edge_sub(1:100:end, 3));

%% Estimate the centroid position 
num_trial = 5; % Number of rows
num_pulse_list = 5 : -1 : 1;
num_p = numel(num_pulse_list); % Number of columns 
x_vec = (0 : 1 : (num_p - 1)) * pulse_spacing_um / pixel_size;
x_vec = round(x_vec - x_vec(ceil(num_trial/2)));
y_vec = (0 : 1 : (num_trial - 1)) * pulse_spacing_um / pixel_size;
y_vec = round(y_vec - y_vec(ceil(num_trial/2)));
init_pos = ceil(im_size(1:2) / 2);
x_vec = x_vec + init_pos(1);
y_vec = y_vec + init_pos(2);

[pulse_x, pulse_y] = meshgrid(x_vec, y_vec);
pulse_dist_to_surf = pulse_x .* n_vec(2) + pulse_y .* n_vec(1) + n_vec(3);

fig_hdl = figure;
ax_hdl = subplot(1,2,1);
imagesc(ax_hdl, plt_z);
ax_hdl.DataAspectRatio = [1,1,1];
hold(ax_hdl, 'on');
scatter(ax_hdl, pulse_x(:), pulse_y(:), 'rx');
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
cbar_hdl = colorbar(ax_hdl);
implay(exp_im_avg{1});
ax_hdl_2 = subplot(1,2,2);
imagesc(ax_hdl_2, bg_mask_zcount);
ax_hdl_2.DataAspectRatio = [1,1,1];
cbar_hdl = colorbar(ax_hdl_2);
%% Find the centroids first 
im_diff_zxy = permute(im_diff, [3, 2, 1]);
im_diff_zxy_avg = mean(im_diff_zxy, 3);

int_peak_x_Q = fun_array_local_maximum(im_diff_zxy_avg, round(pulse_spacing_um * 0.8 / pixel_size));
int_peak_x_idx = find(int_peak_x_Q);
int_peak_x_int = im_diff_zxy_avg(int_peak_x_Q);
int_peak_x_Q = int_peak_x_int >= max(int_peak_x_int) * 0.4; 
int_peak_x_idx = int_peak_x_idx(int_peak_x_Q);
int_peak_sub = fun_ind2sub(size(im_diff_zxy_avg), int_peak_x_idx);
% int_peak_x = max(im_diff_zxy_avg, [], 1);
% int_peak_x_Q = fun_array_local_maximum(int_peak_x, round(pulse_spacing_um * 0.8 / pixel_size));
% int_peak_x_idx = find(int_peak_x_Q);
% int_peak_x_int = int_peak_x(int_peak_x_Q);
% int_peak_x_Q = int_peak_x_int >= max(int_peak_x_int) * 0.75; 
% int_peak_x_idx = int_peak_x_idx(int_peak_x_Q);

% Check 
fig_hdl = figure;
fig_hdl.Position(4) = fig_hdl.Position(4) * 2;
ax_2 = subplot(2,1,1);
vis_slice_idx = round(median(int_peak_sub(:, 1)));
imagesc(ax_2, im_diff(:, :, vis_slice_idx));
ax_2.DataAspectRatio = [1,1,1];
colormap(ax_2, 'gray');
ax_2.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_2.XTick, 'UniformOutput', false);
ax_2.YTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_2.YTick, 'UniformOutput', false);
ax_2.YLabel.String = 'Y (\mum)';
ax_2.XLabel.String = 'X (\mum)';
ax_2.Title.String = sprintf('Slice %d', vis_slice_idx);

ax_1 = subplot(2,1, 2);
imagesc(ax_1, im_diff_zxy_avg);
colormap(ax_1, 'gray');
ax_1.DataAspectRatio = [1,0.25,1];
hold(ax_1, 'on');
ax_1.YLabel.String = 'Z (\mum)';
ax_1.XLabel.String = 'X (\mum)';
ax_1.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_1.XTick, 'UniformOutput', false);
scatter(ax_1, int_peak_sub(:, 2), int_peak_sub(:, 1), 'rx');
ax_1.Title.String = 'Mean projection along Y';

tmp_fig_folder = fullfile(vis_root_folder, 'Site_ablation', site_folder);
tmp_fig_prefix = sprintf('%s_%s', vis_file_prefix, site_folder);
tmp_fig_fp = fullfile(tmp_fig_folder, sprintf('%s_image.png', tmp_fig_prefix));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);

%% Analyze profile one by one 
num_test = numel(int_peak_x_idx);
site_stat = cell(num_test, 1);
for i_peak = 1 : num_test
    tmp_ctr_sub = [ctr_y, int_peak_sub(i_peak, 2), int_peak_sub(i_peak, 1)];
    tmp_im = crop_center_box(im_diff, tmp_ctr_sub, [est_fwhm_y_pxl, est_fwhm_x_pxl, est_fwhm_z_pxl * 3]);
    tmp_fig_folder = fullfile(vis_root_folder, 'Site_ablation', site_folder, num2str(i_peak));
    tmp_fig_prefix = sprintf('%s_%s_t_%d', vis_file_prefix, site_folder, i_peak);
    site_stat{i_peak} = WBIMAnalysisAblationTest.analyze_single_ablation_site(tmp_im, 'pixel_size_um', [0.25, 0.25, 1], ...
        'visQ', true, 'save_fig_Q', true, 'save_folder', tmp_fig_folder, 'save_im_prefix', tmp_fig_prefix);
end

site_stat = cat(1, site_stat{:});
field_list = {'x', 'y', 'z'};
fwhm_stat = [];
fwhm_stat = struct;
for i = 1 : numel(field_list)
    tmp_fn = field_list{i};
    tmp_data = arrayfun(@(x) x.(tmp_fn).FWHM, site_stat);
    fwhm_stat.(tmp_fn) = fun_analysis_get_basic_statistics(tmp_data);
    fwhm_stat.(tmp_fn).data = tmp_data;
end
tmp_data_fp = fullfile(tmp_fig_folder, sprintf('%s_FWHM_stat.mat', tmp_fig_prefix));
save(tmp_data_fp, '-struct', 'fwhm_stat');
fun_save_structure_as_json_file(fwhm_stat, strrep(tmp_data_fp, 'mat', 'json'));
%%
% implay(rescale(max(0, im_diff)))
% Crop image
x_range = round(4 / pixel_size * [-1, 1] + im_size(2) / 2);
y_range = round(75 / pixel_size * [-1, 1] + im_size(1) / 2);
% Remove the first frame. Different laser intensity? 
im_c = im_diff(y_range(1) : y_range(2) - 1, x_range(1) : x_range(2) - 1, :);
    %% Measure the z and x width 
im_c_zyx = permute(im_c, [3, 1, 2]);
im_c_zxy = permute(im_c, [3, 2, 1]);

% Find the center on zx plane by mean projection 
avg_int_zx = mean(im_c_zxy, 3);
% Find the peak 
[max_int, max_ind] = max(avg_int_zx, [], 'all', 'linear');
max_sub = fun_ind2sub(size(avg_int_zx), max_ind);
% Fit gaussian 
x_profile = avg_int_zx(max_sub(1), :);
x_pos_um = (1 : numel(x_profile)) * pixel_size;
fit_para = WBIMAnalysisAblationTest.fit_1d_int_profile_with_gaussian(...
    x_pos_um, x_profile, 'visQ', true, 'vis_x_label', 'X (\mum)');
fig_fp = fullfile(vis_root_folder, 'Single_pulse', sprintf('%s_%s_single_pulse_psfx.png', ...
    vis_file_prefix, 'Site_41'));
fun_print_image_in_several_formats(fit_para.fig_hdl, fig_fp);

z_profile = avg_int_zx(:, max_sub(2));
z_pos_um = 1 : numel(z_profile);
fit_para_z = WBIMAnalysisAblationTest.fit_1d_int_profile_with_gaussian(...
    z_pos_um, z_profile, 'visQ', true, 'vis_x_label', 'Z (\mum)');
fig_fp = fullfile(vis_root_folder, 'Single_pulse', sprintf('%s_%s_single_pulse_psfz.png', ...
    vis_file_prefix, 'Site_41'));
fun_print_image_in_several_formats(fit_para_z.fig_hdl, fig_fp);

% Average around the focal plane: 
z_avg_half_width = ceil(fit_para_z.FWHM/2);
z_avg_range = max_sub(1) + [-1, 1] * z_avg_half_width;
x_avg_range = max_sub(2) + [-1, 1] * ceil(fit_para.FWHM/2/pixel_size);
avg_int_yx = mean(im_c(:, :, z_avg_range(1) : z_avg_range(2)), 3);
y_profile = mean(avg_int_yx(:, x_avg_range(1) : x_avg_range(2)), 2);
fit_para_y = WBIMAnalysisAblationTest.fit_1d_int_profile_with_gaussian(...
    (1 : numel(y_profile)) * pixel_size, y_profile, 'visQ', true, ...
    'vis_x_label', 'Y (\mum)');
fig_fp = fullfile(vis_root_folder, 'Single_pulse', sprintf('%s_%s_single_pulse_psfy.png', ...
    vis_file_prefix, 'Site_41'));
fun_print_image_in_several_formats(fit_para_y.fig_hdl, fig_fp);

fig_hdl = figure;
ax_1 = subplot(1,2,1);
imagesc(avg_int_zx);
ax_1.DataAspectRatio = [1,0.25,1];
ax_1.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_1.XTick, 'UniformOutput', false);
ax_1.XLabel.String = 'X (\mum)';
ax_1.YLabel.String = 'Z (\mum)';

ax_2 = subplot(1,2,2);
imagesc(avg_int_yx);
ax_2.DataAspectRatio = [1,1,1];
ax_2.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_2.XTick, 'UniformOutput', false);
ax_2.YTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_2.YTick, 'UniformOutput', false);
ax_2.XLabel.String = 'X (\mum)';
ax_2.YLabel.String = 'Y (\mum)';
fig_fp = fullfile(vis_root_folder, 'Single_pulse', sprintf('%s_%s_single_pulse_mean_proj.png', ...
    vis_file_prefix, 'Site_41'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

tmp = WBIMAnalysisAblationTest.analyze_single_ablation_site(im_c, 'pixel_size_um', [0.25, 0.25, 1], ...
    'visQ', true);


% fig_hdl = figure;
% ax_hdl = subplot(1,3,1);
% imagesc(im_c(:, :, 52));
% ax_hdl.DataAspectRatio = [1,1,1];
% ax_2 = subplot(1,3,2);
% imagesc(avg_int_yx);
% ax_2.DataAspectRatio = [1,1,1];
% ax_3 = subplot(1,3,3);
% plot(ax_3, y_profile);


%%
vis_frame_idx = 271 : 275;
vis_stack = exp_im{2}(:, :, vis_frame_idx);
vis_stack = vis_stack(y_range(1) : y_range(2) - 1, x_range(1) : x_range(2) - 1, :);
num_vis = numel(vis_frame_idx);
avg_im = uint16(median(vis_stack, 3));

fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
for i_vis = 1 : num_vis
    ax_hdl = subplot(1, num_vis + 1, i_vis);
    tmp_vis_im = vis_stack(:, :, i_vis) - avg_im;
    imagesc(tmp_vis_im);
    ax_hdl.DataAspectRatio = [1,1,1];
%     ax_hdl.Title.String = sprintf('%d s', round(vis_frame_idx(i_vis) / frame_rate_Hz + 9.7));
%     ax_hdl.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_hdl.XTick, 'UniformOutput', false);
%     ax_hdl.YTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_hdl.YTick, 'UniformOutput', false);
%     colormap(ax_hdl, 'gray');
%     if i_vis == 1
%         ax_hdl.XLabel.String = 'X (\mum)';
%         ax_hdl.YLabel.String = 'Y (\mum)';
%     end
end
ax_hdl = subplot(1, num_vis + 1, num_vis + 1);
imagesc(avg_im);
ax_hdl.DataAspectRatio = [1,1,1];
%% Motion correction 
tmp_ref = avg_im;
tmp_reg = vis_stack(:, :, 5);
tmp_ref_mask = tmp_ref >= 0;
tmp_reg_mask = tmp_reg >= 0;
fft_search_range = [10, 10];
[disp_vec, tmp2] = MaskedTranslationRegistration(tmp_ref, tmp_reg, ...
    tmp_ref_mask, tmp_reg_mask, fft_search_range);
%% Row shift detection and correction 
tmp_sec = 2;
local_im = vis_stack(:, :, tmp_sec);
[tmp_shift, im_c] = fun_estimate_2p_im_line_shift(local_im, 20);

%%









