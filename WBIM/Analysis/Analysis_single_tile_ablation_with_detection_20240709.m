clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'Ablation_test';
% site_idx = 5;
% exp_name = sprintf('WBIMSkl_20240305002_%s_%d', '20240710', site_idx);
exp_name = sprintf('WBIM20240305002_%s', '20240710');
scan_tiles_cell = DataManager.load_tile_in_experiment(exp_group, exp_name, WBIMMicroscopeMode.Scan);
explore_tiles_cell = DataManager.load_tile_in_experiment(exp_group, exp_name, WBIMMicroscopeMode.Explore);
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'Damage_range');
%%
scan_tile_info = scan_tiles_cell(1:end-1);
scan_tile = cellfun(@(x) x.load_tile(), scan_tile_info, 'UniformOutput', false);
scan_tile_RBG = cellfun(@(x) fun_merge_image_stacks(x, 'method', 'GraySkull', 'stretch_contrast_Q', true), ...
    scan_tile, 'UniformOutput', false);
%%
vis_layer = 5;
implay(scan_tile_RBG{vis_layer});
%% 
l1_vis_sec = 199;
vis_ly_idx = [1, 2, 3, 4, 5, 6, 7, 8];
num_vis = numel(vis_ly_idx);
layer_overlap_um = 25;
fig_hdl = figure;
for i = 1 : num_vis
    ax = nexttile();
    vis_sec = l1_vis_sec - (vis_ly_idx(i) - vis_ly_idx(1)) * layer_overlap_um;
    
    image(ax, scan_tile_RBG{vis_ly_idx(i)}(:, :, :, vis_sec));
    ax.DataAspectRatio = [1,1,1];
    ax.Title.String = sprintf('Layer %d, z = %d \\mum', vis_ly_idx(i), vis_sec);
end
%%
all_tiles = cat(1, scan_tiles_cell{:}, explore_tiles_cell{:});
% Sort by time: 
tile_acq_time = cat(1, all_tiles.t_done);
[tile_acq_time, t_idx] = sortrows(tile_acq_time);
all_tiles = all_tiles(t_idx);
%%
scan_tile_idx = 6;
scan_tile_t_idx = find(t_idx == scan_tile_idx);
scan_tile_tp1_idx = find(t_idx == (scan_tile_idx + 1));
exp_tile_after_scan_idx = (scan_tile_t_idx + 1) : (scan_tile_tp1_idx - 1);
exp_tile_after_scan_info = all_tiles(exp_tile_after_scan_idx);
exp_tile_after_scan = arrayfun(@(x) x.load_tile(), exp_tile_after_scan_info, 'UniformOutput', false);
exp_tile_after_scan_RGB = cellfun(@(x) fun_merge_image_stacks(x, 'method', 'GraySkull', 'stretch_contrast_Q', true), ...
    exp_tile_after_scan, 'UniformOutput', false);
%%
abl_sec_idx = 1 : 25;
s_rgb_abl_mip = squeeze(max(scan_tile_RBG{scan_tile_idx}(:, :, :, abl_sec_idx), [], 4));
e_rgb_abl_mip = squeeze(max(exp_tile_after_scan_RGB{1}(:, :, :, abl_sec_idx), [], 4));
s_skl_abl_mip = max(scan_tile{scan_tile_idx}{2}(:, :, abl_sec_idx), [], 3);
s_vsl_abl_mip = max(scan_tile{scan_tile_idx}{1}(:, :, abl_sec_idx), [], 3);
e_skl_abl_mip = max(exp_tile_after_scan{1}{2}(:, :, abl_sec_idx), [], 3);
e_vsl_abl_mip = max(exp_tile_after_scan{1}{1}(:, :, abl_sec_idx), [], 3);
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 2;
ax_1 = nexttile();
image(ax_1, s_rgb_abl_mip);
ax_1.DataAspectRatio = [1,1,1];
ax_1.Title.String = sprintf('Layer %d Scan, z [%d, %d] \\mum', scan_tile_idx, abl_sec_idx([1,end]));
ax_2 = nexttile();
image(ax_2, e_rgb_abl_mip);
ax_2.DataAspectRatio = [1,1,1];
ax_2.Title.String = sprintf('Layer %d Explore 1, z [%d, %d] \\mum', scan_tile_idx, abl_sec_idx([1,end]));
ax_3 = nexttile();
imshowpair(s_skl_abl_mip, e_skl_abl_mip);
ax_3.Title.String = sprintf('Skull overlay');
ax_4 = nexttile();
imshowpair(s_vsl_abl_mip, e_vsl_abl_mip);
ax_4.Title.String = sprintf('Vessel overlay');
fig_fp = fullfile(vis_folder, sprintf('%s_%s_scan_6_vs_exp_7_1.png', exp_group, exp_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% 
s_rgb_zxy = permute(scan_tile_RBG{scan_tile_idx}, [4, 2, 3, 1]);
vis_sec_idx = 150 : 100 : 550;

fig_hdl = figure;
fig_hdl.Position(4) = fig_hdl.Position(4) * 4;
t = tiledlayout(numel(vis_sec_idx), 1);
for i = vis_sec_idx
    ax_1 = nexttile();
    image(ax_1, s_rgb_zxy(:, :, :, i));
    ax_1.DataAspectRatio = [1,1,1];
    ax_1.Title.String = sprintf('y = %d \\mum', i);
    ax_1.XLabel.String = 'X (\mum)';
    ax_1.YLabel.String = 'Z (\mum)';
end
fig_fp = fullfile(vis_folder, sprintf('%s_%s_scan_6_zx_selected_y_sec.png', exp_group, exp_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

