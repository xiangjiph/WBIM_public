clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'Ablation_test';
exp_name = 'WBIMGel_20240306';
data_root_folder = fullfile(DataManager.fp_layer(exp_group, exp_name, 1), ...
    'Explore', 'Continuous_ablation');
si_im_cvrt = @(x) 2 * uint16(x);

vis_root_folder = fullfile(DataManager.fp_visualization_folder(exp_group, exp_name), ...
    'Continuous_ablation');
vis_file_prefix = sprintf('%s_%s_cont_abl', exp_group, exp_name);
%% Iterative imaging and ablation
site_idx = 22;
site_folder = sprintf('Site_%d', site_idx);
exp_folder = fullfile(data_root_folder, site_folder);
exp_files = dir(fullfile(exp_folder, '*.tif'));

exp_fig_folder = fullfile(vis_root_folder, site_folder);
exp_fig_prefix = sprintf('%s_cont_abl_%s', vis_file_prefix, site_folder);

file_name_list = {exp_files.name};
file_count_list = cellfun(@(x) str2double(x(6:10)), file_name_list);
assert(issorted(file_count_list, 'ascend'));
num_files = numel(file_count_list);
% Load all files 
exp_im = arrayfun(@(x) si_im_cvrt(DataManager.load_single_tiff(fullfile(exp_folder, ...
    x.name))), exp_files, 'UniformOutput', false);

pixel_size = 1;
frame_avg = 5;

% est_fwhm_z_um = 5;
% est_fwhm_y_um = 75;
% est_fwhm_x_um = 4;
% 
% est_fwhm_y_pxl = round(est_fwhm_y_um / pixel_size);
% est_fwhm_x_pxl = round(est_fwhm_x_um / pixel_size);
% est_fwhm_z_pxl = round(est_fwhm_z_um);
% ctr_y = round(im_size(1) / 2);

registration_Q = true;
%%
% Smooth the image using median filter along the y axis
% exp_im_sm = cellfun(@(x) medfilt3(x, [3, 3, 1]), exp_im, 'UniformOutput', false);
exp_im_avg = cellfun(@(x) fun_image_stack_average_frame(x, frame_avg, @mean), ...
    exp_im, 'UniformOutput', false);
im_size = size(exp_im_avg{1});
% Move the tile after ablation up? 
%% Registration - for potential global deformation 
% Global translation is relatively small (~ 1 um in y and z). Neglect them
% at the moment.
if registration_Q
    mfft_th = 7e4;
    fft_search_range = [5, 5, 10];
    
    ref_im = exp_im_avg{1};
    for i = 2 : numel(exp_im_avg)
        tmp_reg = exp_im_avg{i};
        tmp_ref_mask = ref_im >= mfft_th;
        tmp_reg_mask = tmp_reg >= mfft_th;
        [tmp_vec, tmp2] = MaskedTranslationRegistration(ref_im, tmp_reg, ...
            tmp_ref_mask, tmp_reg_mask, fft_search_range);
        if any(tmp_vec)
            tmp_reg_mv = imtranslate(tmp_reg, tmp_vec);
            tmp_empty_mask = (tmp_reg_mv == 0);
            tmp_reg_mv(tmp_empty_mask) = tmp_reg(tmp_empty_mask);
            
            % Check the zy profile
            fig_hdl = figure;
            ax_1 = subplot(2,1,1);
            imshowpair(squeeze(ref_im(:, 300, :)).', squeeze(tmp_reg(:, 300, :)).');
            ax_1.Title.String = sprintf('Stack %d Displacement vector [%d, %d, %d]',...
                i, tmp_vec);
            ax_2 = subplot(2,1,2);
            imshowpair(squeeze(ref_im(:, 300, :)).', squeeze(tmp_reg_mv(:, 300, :)).');
            exp_im_avg{i} = tmp_reg_mv;
            fprintf('Stack %d Displacement vector [%d, %d, %d]\n', i, tmp_vec);
        end
    end
end
% im_diff = exp_im_avg{1} - exp_im_avg{2};
% im_size = size(im_diff);
exp_im_avg_sm = cellfun(@(x) medfilt3(x, [3, 3, 3]), exp_im_avg, 'UniformOutput', false);
%%
exp_im_avg_diff = cellfun(@(x1, x2) x1 - x2, exp_im_avg_sm(1:end-1), exp_im_avg_sm(2:end), 'UniformOutput', false);
% Compute the baseline surface height 
im_baseline = exp_im_avg_sm{1};
otsu_th = graythresh(im_baseline) * double(intmax(class(im_baseline)));
tmp_mask = im_baseline > otsu_th;
surface_dist_map = sum(~tmp_mask, 3);

row_shift_col = [1 : 10, (-9 : 0) + im_size(2)];
surface_dist_map(:, row_shift_col) = nan;

abl_idx = 1;
tmp_diff = exp_im_avg_diff{abl_idx};
% diff_mask_th = otsu_th;
diff_mask_th = 1e4;
tmp_diff_mask = tmp_diff > diff_mask_th;
tmp_diff_mask_mip = any(tmp_diff_mask, 3);

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
im_hdl = contour(ax_hdl, surface_dist_map, 'ShowText', 'on');
ax_hdl.YDir = 'reverse';
hold(ax_hdl, 'on');
im_2 = imagesc(ax_hdl, tmp_diff_mask_mip * -1, 'AlphaData', tmp_diff_mask_mip * 0.75);
ax_hdl.DataAspectRatio = [1,1,1];
c_bar = colorbar(ax_hdl);
cmap = colormap(ax_hdl,'winter'); 
cmap(1, :) = [1, 0, 0];
ax_hdl.Colormap = cmap;
ax_hdl.CLim(1) = 0;
c_bar.Limits(1) = 1;
c_bar.Label.String = 'Distance to the surface (\mum)';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
rectangle(ax_hdl, 'Position', [100, 50, 500, 600], 'LineStyle', '--', 'EdgeColor', 'k');
ax_hdl.Title.String = sprintf('Intensity density after ablation %d', abl_idx);

tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_surface_damage_vis_abl_%d.png', exp_fig_prefix, abl_idx));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);

%% Look at the max projection of the image difference 
vis_x_idx = 150 : 550;
vis_mip_zyx = squeeze(max(im_baseline(:, vis_x_idx, :), [], 2));
vis_diff_mip_zyx = squeeze(max(exp_im_avg_diff{abl_idx}(:, vis_x_idx, :), [], 2));

vis_im_rbg = repelem(vis_mip_zyx, 1, 1, 3);
vis_im_rbg(:, :, 1) = max(vis_im_rbg(:, :, 1), vis_diff_mip_zyx * 3);
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
image(ax_hdl, permute(vis_im_rbg, [2, 1, 3]));
ax_hdl.DataAspectRatio = [1,1,1];
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Y (\mum)';
ax_hdl.YLabel.String = 'Z (\mum)';
ax_hdl.Title.String = sprintf('MIP X in [%d, %d] \\mum', vis_x_idx([1, end]));
% imshowpair(vis_mip_zyx.', vis_diff_mip_zyx.');
tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_abl_%d_diff_vis_x_%d_%d.png',...
    exp_fig_prefix, abl_idx, vis_x_idx([1, end])));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);

%% Difference 
vis_mip_yx = max(im_baseline, [], 3);
vis_diff_mip_yx = max(exp_im_avg_diff{abl_idx}, [], 3);

vis_im_rbg = repelem(vis_mip_yx, 1, 1, 3);
vis_im_rbg(:, :, 1) = max(vis_im_rbg(:, :, 1), vis_diff_mip_yx * 1.5);
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
image(ax_hdl, vis_im_rbg);
ax_hdl.DataAspectRatio = [1,1,1];
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.Title.String = sprintf('MIP Z');
tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_abl_%d_diff_vis_z.png',...
    exp_fig_prefix, abl_idx));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
%% After ablation zxy
vis_im_zxy = permute(exp_im_avg_sm{end}, [3, 2, 1]);
vis_sec_list = 200 : 100 : 600;
num_vis = numel(vis_sec_list);

fig_hdl = figure;
fig_hdl.Position(4) = fig_hdl.Position(4) * 2;
% ax_hdl = axes(fig_hdl);
for i = 1 : num_vis
    vis_sec = vis_sec_list(i);
    ax_hdl = subplot(num_vis, 1, i);
    imagesc(ax_hdl, vis_im_zxy(:, :, vis_sec));
    colormap(ax_hdl, 'gray');
    ax_hdl.DataAspectRatio = [1,1,1];
    hold(ax_hdl, 'on');
    rectangle(ax_hdl, 'Position', [100, 0, 500, 50], 'LineStyle', '--', 'EdgeColor', 'r');
    ax_hdl.XLabel.String = 'Z (\mum)';
    ax_hdl.YLabel.String = 'X (\mum)';
    ax_hdl.Title.String = sprintf('Y = %d \\mum', vis_sec);
end

tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_abl_%d_aft_vis_zxy_secs.png',...
    exp_fig_prefix, abl_idx));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);

%%
exp_im_avg_zyx = cellfun(@(x) permute(x, [3, 1, 2]), exp_im_avg_sm, 'UniformOutput', false);
vis_sec_list = [200, 300, 400, 500];
num_sec = numel(vis_sec_list);
compare_list = 2 : numel(exp_im_avg_zyx);
num_vis_col = numel(compare_list);
fig_hdl = figure;
ref_stack = exp_im_avg_zyx{1};
for i_stack = 1 : num_vis_col
    tmp_stack = exp_im_avg_zyx{compare_list(i_stack)};
    for i_sec = 1 : num_sec
        ax = subplot(num_sec, num_vis_col, sub2ind([num_vis_col, num_sec], i_stack, i_sec));
        tmp_sec = vis_sec_list(i_sec);
        imshowpair(ref_stack(:, :, tmp_sec), tmp_stack(:, :, tmp_sec));        
        ax.Title.String = sprintf('X = %d \\mum', tmp_sec);
        if i_sec == num_sec
           ax.XLabel.String = 'Y (\mum)';
           ax.XAxis.Visible = 'on';
           ax.YAxis.Visible = 'on';
           ax.YLabel.String = 'Z (\mum)';
        end
        hold(ax, 'on');
%         line(ax, [50, 650], [50, 50], 'Color', 'r', 'LineStyle', '--');
        rectangle(ax, 'Position', [100, 0, 500, 50], 'LineStyle', '--', 'EdgeColor', 'r');
    end    
end
tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_abl_%d_aft_vis_zyx_compare_secs.png',...
    exp_fig_prefix, abl_idx));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);


%% Visualize damage
compare_list = 1 : 4;
num_vis = numel(compare_list);
fig_hdl = figure;
for i_vis = compare_list
    ax_hdl = subplot(num_vis, 1, i_vis);
    tmp_diff = exp_im_avg_zyx{i_vis} - exp_im_avg_zyx{i_vis+1};
    tmp_diff_max = max(tmp_diff, [], 3);
    
    imagesc(ax_hdl, tmp_diff_max);
    cbar_hdl = colorbar(ax_hdl);
    
    ax_hdl.CLim(1) = 1e4;
    ax_hdl.DataAspectRatio = [1,1,1];
    if i_vis == 1
        ax_hdl.Title.String = sprintf('MIP of intensity difference with the previous tack');
    end
    if i_vis == num_vis
        ax_hdl.XLabel.String = 'Y (\mum)';
        ax_hdl.YLabel.String = 'Z (\mum)';
    end
end
%%
volumeViewer(exp_im_avg_zyx{1} - exp_im_avg_zyx{4})



%% Iterative plane ablation visualization: depth vs ablation plane
im_size = size(exp_im_avg{1});
vis_plane_idx = 1 + (0 : (num_files -1)) * 5;
vis_dz_um = vis_plane_idx - 1;
abl_file_idx = 1 : num_files;
num_vis_stack = numel(abl_file_idx);
num_vis_plane = numel(vis_plane_idx);
% Bounding box
tmp_bbox_ctr = round(im_size(1:2) / 2);
tmp_bbox_ll_um = [100, 50];
tmp_bbox_ll = tmp_bbox_ll_um / pixel_size;
tmp_bbox_mmll = [tmp_bbox_ctr - tmp_bbox_ll/2, tmp_bbox_ll];

tic
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 3;
for i_stack = 1 : num_vis_stack
    tmp_stack = exp_im_avg{abl_file_idx(i_stack)};
    for i_plane = 1 : num_vis_stack
        tmp_plane_idx = vis_plane_idx(i_plane);
        tmp_im = tmp_stack(:, :, tmp_plane_idx);
        tmp_ind = sub2ind([num_vis_stack, num_vis_plane], i_stack, i_plane);
        tmp_ax = subplot(num_vis_plane, num_vis_stack, tmp_ind);        
%         imshow(tmp_im);
        imagesc(tmp_ax, tmp_im, 'CDataMapping', 'direct');
        tmp_ax.DataAspectRatio = [1,1,1];
        colormap(tmp_ax, 'gray');
        tmp_ax.XTick = [];
        tmp_ax.YTick = [];
        if i_plane == (i_stack - 1)
            hold(tmp_ax, 'on');
            rectangle(tmp_ax, 'Position', tmp_bbox_mmll, 'LineStyle', '--', 'EdgeColor', 'r');
        end        
        if i_stack == 1 && i_plane == 1
            tmp_ax.Title.String = 'Before';
        elseif i_plane == 1
            tmp_ax.Title.String = sprintf('Ablation @ %d \\mum', vis_dz_um(i_stack - 1));
        end
        if i_stack == 1
            tmp_ax.YLabel.String = sprintf('Z = %d \\mum', 5 * (i_plane - 1));
        end
    end
end
toc
exp_fig_folder = fullfile(vis_root_folder, 'Site_ablation', site_folder);
exp_fig_prefix = sprintf('%s_%s', vis_file_prefix, site_folder);
tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_image.png', exp_fig_prefix));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
%% Volume ablation
exp_im_avg_zxy = cellfun(@(x) permute(x, [3, 2, 1]), exp_im_avg, 'UniformOutput', false);
abl_param = WBIMSPAblation.compute_ablation_parameters(1.5, 0.75, false, 5);
fast_buffer_um = abl_param.fast_axis_speed_um_s * WBIMConfig.ABLATION_ZABER_PRETRIGGER_TIME_UNCERTAINTY_ms / ...
    1e3;

vis_sec = 150;
selected_x_idx = 50 : 650;
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
ax_1 = subplot(2,1,1);
imagesc(ax_1, exp_im_avg_zxy{1}(:, selected_x_idx, vis_sec));
ax_1.DataAspectRatio = [1, 1, 1];
ax_1.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_1.XTick, 'UniformOutput', false);
ax_1.Title.String = 'Before';
ax_1.XLabel.String = 'X (\mum)';
ax_1.YLabel.String = 'Z (\mum)';
colormap(ax_1, 'gray');

ax_2 = subplot(2,1,2);
imagesc(ax_2, exp_im_avg_zxy{2}(:, selected_x_idx, vis_sec));
colormap(ax_2, 'gray');
hold(ax_2, 'on');
abl_x_range = ((im_size(2) / 2) + [-258, 258] / pixel_size) - selected_x_idx(1);
for i = 0 : 5 : 50 
    line(ax_2, abl_x_range, [i, i], 'Color', 'r', 'LIneStyle', '--');    
end
ax_2.YLim(1) = 0.5;
ax_2.DataAspectRatio = [1,1,1];
ax_2.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_2.XTick, 'UniformOutput', false);
ax_2.Title.String = '(1.5, 0.75, 5 J/cm^2, 5 \mum)';
ax_2.XLabel.String = 'X (\mum)';
ax_2.YLabel.String = 'Z (\mum)';

exp_fig_folder = fullfile(vis_root_folder, site_folder);
exp_fig_prefix = sprintf('%s_%s', vis_file_prefix, site_folder);
tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_zx_sec_%d.png', exp_fig_prefix, vis_sec));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
%%
tmp_bbox_ctr = round(im_size(1:2) / 2);
tmp_bbox_ll_um = [500, 500];
tmp_bbox_ll = tmp_bbox_ll_um / pixel_size;
tmp_bbox_mmll = [tmp_bbox_ctr - tmp_bbox_ll/2, tmp_bbox_ll];

vis_sec_idx = [0, 25, 50, 55, 60] + 1;
num_sec = numel(vis_sec_idx);
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 4;
for i = 1 : num_sec
    ax_1 = subplot(2, num_sec, sub2ind([num_sec, 2], i, 1));
    ax_2 = subplot(2, num_sec, sub2ind([num_sec, 2], i, 2));
    imagesc(ax_1, exp_im_avg{1}(:, :, vis_sec_idx(i)));
    imagesc(ax_2, exp_im_avg{end}(:, :, vis_sec_idx(i)));
    colormap(ax_1, 'gray');
    colormap(ax_2, 'gray');
    [ax_1.DataAspectRatio, ax_2.DataAspectRatio] = deal([1,1,1]);
    if vis_sec_idx(i) <= 51
        hold(ax_2, 'on');
        rectangle(ax_2, 'Position', tmp_bbox_mmll, 'LineStyle', '--', 'EdgeColor', 'r');
    end
    if i == 1
        ax_1.YLabel.String = 'Before';
        ax_2.YLabel.String = 'After';
        ax_1.XTick = [];
        [ax_1.YTick, ax_2.YTick] = deal([]);
        ax_2.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_2.XTick, 'UniformOutput', false);
        ax_2.XLabel.String = 'X (\mum)';
    else
        [ax_1.XTick, ax_2.XTick] = deal([]);
        [ax_1.YTick, ax_2.YTick] = deal([]);
    end
    ax_1.Title.String = sprintf('%d \\mum', vis_sec_idx(i) - 1);
end

exp_fig_folder = fullfile(vis_root_folder, site_folder);
exp_fig_prefix = sprintf('%s_%s', vis_file_prefix, site_folder);
tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_yx_5z.png', exp_fig_prefix));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);