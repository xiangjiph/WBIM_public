clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'Ablation_test';
exp_name = 'WBIMGel_20240306';
data_root_folder = fullfile(DataManager.fp_layer(exp_group, exp_name, 1), ...
    'Explore');
si_im_cvrt = @(x) 2 * uint16(x);

vis_root_folder = DataManager.fp_visualization_folder(exp_group, exp_name);
vis_file_prefix = sprintf('%s_%s', exp_group, exp_name);
%% Iterative imaging and ablation
site_idx = 48;
site_folder = sprintf('Site_%d', site_idx);
exp_folder = fullfile(data_root_folder, site_folder);
exp_files = dir(fullfile(exp_folder, '*.tif'));
file_name_list = {exp_files.name};
file_count_list = cellfun(@(x) str2double(x(6:10)), file_name_list);
assert(issorted(file_count_list, 'ascend'));
num_files = numel(file_count_list);
% Load all files 
exp_im = arrayfun(@(x) si_im_cvrt(DataManager.load_single_tiff(fullfile(exp_folder, ...
    x.name))), exp_files, 'UniformOutput', false);

pixel_size = 0.25;
frame_avg = 5;

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
% Global translation is relatively small (~ 1 um in y and z). Neglect them
% at the moment.
if registration_Q
    tmp_ref = exp_im_avg{1};
    tmp_reg = exp_im_avg{3};
    tmp_ref_mask = tmp_ref >= 1e4;
    tmp_reg_mask = tmp_reg >= 1e4;
    fft_search_range = [10, 10, 10];
    [tmp1, tmp2] = MaskedTranslationRegistration(tmp_ref, tmp_reg, ...
        tmp_ref_mask, tmp_reg_mask, fft_search_range);
    if any(tmp1)
        tmp_reg_mv = imtranslate(tmp_reg, tmp1);
        tmp_empty_mask = (tmp_reg_mv == 0);
        tmp_reg_mv(tmp_empty_mask) = tmp_reg(tmp_empty_mask);
        
        % Check the zy profile
        fig_hdl = figure;
        ax_1 = subplot(1,2,1);
        imshowpair(squeeze(tmp_ref(:, 300, :)), squeeze(tmp_reg(:, 300, :)));
        ax_2 = subplot(1,2,2);
        imshowpair(squeeze(tmp_ref(:, 300, :)), squeeze(tmp_reg_mv(:, 300, :)));
%         exp_im_avg{2} = tmp_reg_mv;
    end
end
% im_diff = exp_im_avg{1} - exp_im_avg{2};
% im_size = size(im_diff);
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
tmp_fig_folder = fullfile(vis_root_folder, 'Site_ablation', site_folder);
tmp_fig_prefix = sprintf('%s_%s', vis_file_prefix, site_folder);
tmp_fig_fp = fullfile(tmp_fig_folder, sprintf('%s_image.png', tmp_fig_prefix));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
%% Volume ablation
exp_im_avg_zxy = cellfun(@(x) permute(x, [3, 2, 1]), exp_im_avg, 'UniformOutput', false);

vis_sec = 329;
selected_x_idx = 50 : 650;
fig_hdl = figure;
ax_1 = subplot(2,1,1);
imagesc(ax_1, exp_im_avg_zxy{1}(:, selected_x_idx, vis_sec));
ax_1.DataAspectRatio = [1, 0.25, 1];
ax_1.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_1.XTick, 'UniformOutput', false);
ax_1.Title.String = 'Before';
ax_1.XLabel.String = 'X (\mum)';
ax_1.YLabel.String = 'Z (\mum)';
colormap(ax_1, 'gray');

ax_2 = subplot(2,1,2);
imagesc(ax_2, exp_im_avg_zxy{2}(:, selected_x_idx, vis_sec));
colormap(ax_2, 'gray');
hold(ax_2, 'on');
abl_x_range = ((im_size(2) / 2) + [-50, 50] / pixel_size) - selected_x_idx(1);
for i = 0.5 : 5 : 45 
    line(ax_2, abl_x_range, [i, i], 'Color', 'r', 'LIneStyle', '--');    
end
ax_2.YLim(1) = 0.5;
ax_2.DataAspectRatio = [1,0.25,1];
ax_2.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_2.XTick, 'UniformOutput', false);
ax_2.Title.String = '(1, 0.75, 5 J/cm^2, 5 \mum)';
ax_2.XLabel.String = 'X (\mum)';
ax_2.YLabel.String = 'Z (\mum)';

tmp_fig_folder = fullfile(vis_root_folder, 'Volume_ablation', site_folder);
tmp_fig_prefix = sprintf('%s_%s', vis_file_prefix, site_folder);
tmp_fig_fp = fullfile(tmp_fig_folder, sprintf('%s_zx_sec_329.png', tmp_fig_prefix));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
%%
tmp_bbox_ctr = round(im_size(1:2) / 2);
tmp_bbox_ll_um = [100, 50];
tmp_bbox_ll = tmp_bbox_ll_um / pixel_size;
tmp_bbox_mmll = [tmp_bbox_ctr - tmp_bbox_ll/2, tmp_bbox_ll];

vis_sec_idx = (0 : 5 : 50) + 1;
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
    if i <= 10
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

tmp_fig_folder = fullfile(vis_root_folder, 'Volume_ablation', site_folder);
tmp_fig_prefix = sprintf('%s_%s', vis_file_prefix, site_folder);
tmp_fig_fp = fullfile(tmp_fig_folder, sprintf('%s_yx_0_5_50.png', tmp_fig_prefix));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);