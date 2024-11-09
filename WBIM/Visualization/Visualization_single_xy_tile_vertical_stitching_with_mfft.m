% Stitching with static lens deformation correction 
clc;clear;close all;
DataManager = WBIMFileManager;
%%
exp_group = 'TestSample';
exp_name = 'WBIM20230425001_20230515';
% DataManager.download_tile_info_in_experiment(exp_group, exp_name);
tile_str = DataManager.load_tile_in_experiment(exp_group, exp_name);
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'Stitched');
%%
stitch_voxel_size_um = [1,1,1];

acq_para_fp = DataManager.fp_acq_control_parameter(exp_group, exp_name, 0);
acq_para = DataManager.load_data(acq_para_fp);
im_grids = acq_para.imaging.grid;
roi_detection_channel = acq_para.parameters.imaging.active_channel;

exp_z0_offset_um = im_grids(WBIMMicroscopeMode.Scan).piezo_z_list_um(1) - ...
    im_grids(WBIMMicroscopeMode.Explore).piezo_z_list_um(1);
exp_dz_um = tiles_in_layer(1).pixel_size_um(3);
num_layers_above_scan_plane = round(exp_z0_offset_um / exp_dz_um);
%% Visualize SMIPs
vis_layer_idx_list = [1,2,3,4,5,6];
num_vis_layer = numel(vis_layer_idx_list);
merged_smip_list = cell(num_vis_layer, 1);
for vis_l = 1 : num_vis_layer
    tmp_tiles = WBIMTileManager.select_duplicated_tiles(tile_str.Scan{vis_layer_idx_list(vis_l)}, 'oldest');
    [tmp_smip_cell, ~] = WBIMTileManager.get_stitched_step_mip(tmp_tiles, roi_detection_channel, 'nanTo0Q', true);
    merged_smip_list{vis_l} = fun_merge_image_stacks(tmp_smip_cell([1,2,4]), 'method', 'GraySkull', 'stretch_contrast_Q', true);
end

implay(merged_smip_list{1});
%% Single tile vertical stitching 
% Find candidate tile:
layer_idx = 2;
tiles_in_layer = WBIMTileManager.select_duplicated_tiles(tile_str.Scan{layer_idx}, 'oldest');
% Get step mip
[smip_cell, local_tile_info] = WBIMTileManager.get_stitched_step_mip(...
    tiles_in_layer, roi_detection_channel, 'nanTo0Q', true);

smip_merge = fun_merge_image_stacks(smip_cell([1,2,4]), 'method', 'GraySkull', 'stretch_contrast_Q', true);
%%
tv = WBIMTileViewer(tiles_in_layer);
tv.update_tile_grid();
tv.init_step_mip(smip_cell{2}, local_tile_info.step_mip_pixel_yxz_um);
% tv.init_step_mip_mask(
tv.ui_init_ctrl();
%%
vis_idx = 177;
vis_tile = tiles_in_layer(vis_idx);
vis_stack_c = vis_tile.load_tile();
vis_stack_c_um = cellfun(@(x) imresize3(fun_stretch_contrast(medfilt3(x)), vis_tile.stack_size_um), ...
    vis_stack_c, 'UniformOutput', false);
vis_stack_yxz = fun_merge_image_stacks(vis_stack_c_um, 'method', 'VesselBone');
vis_stack_zyx = permute(vis_stack_yxz, [4, 1, 3, 2]);
vis_stack_zxy = permute(vis_stack_yxz, [4, 2, 3, 1]);

implay(vis_stack_zxy)
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1.5];
ax1 = subplot(1,3,1);
image(ax1, squeeze(max(vis_stack_yxz, [], 4)));
hold(ax1, 'on');
xline(ax1, 55, 'w-.');
xline(ax1, 35, 'w-.');

yline(ax1, 300, 'w-.');
yline(ax1, 320, 'w-.');

ax1.DataAspectRatio = [1,1,1];
ax1.Title.String = 'Z MIP';
ax1.XLabel.String = 'X (\mum)';
ax1.YLabel.String = 'Y (\mum)';
ax2 = subplot(1,3,2);
image(ax2, permute(squeeze(max(vis_stack_zyx(:, :, :, 35:55), [], 4)), [2, 1, 3]));
ax2.DataAspectRatio = [1,1,1];
ax2.Title.String = 'X 35 - 55 \mum MIP ';
ax2.XLabel.String = 'Z (\mum)';
ax2.YLabel.String = 'Y (\mum)';
ax2.XDir = 'reverse';
ax3 = subplot(1,3,3);
image(ax3, squeeze(max(vis_stack_zxy(:, :, :, 300:320), [], 4)));
ax3.DataAspectRatio = [1,1,1];
ax3.Title.String = 'Y 300 - 320 \mum MIP ';
ax3.XLabel.String = 'X (\mum)';
ax3.YLabel.String = 'Z (\mum)';

fig_fp = fullfile(vis_folder, sprintf('%s_%s_layer_%d_grid_%d_mip.png', exp_group, ...
    exp_name, layer_idx, vis_tile.grid_ind));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
% Get all tiles at the same xy position: 
vis_grid_sub = vis_tile.grid_sub;
vert_tiles = {};
for i = 1 : numel(tile_str.Scan)
    tmp_ts = WBIMTileManager.select_duplicated_tiles(tile_str.Scan{i}, 'latest');
    tmp_subs = cat(1, tmp_ts.grid_sub);
    tmp_Q = all(tmp_subs == vis_grid_sub, 2);
    if any(tmp_Q)
        tmp_t = tmp_ts(tmp_Q);
        assert(numel(tmp_t) == 1);
        vert_tiles{end+1} = tmp_t;
    end
end
vert_tiles = cat(1, vert_tiles{:});
num_vt = numel(vert_tiles);
ref_z_um = vert_tiles(1).sample_xyz_um_done(3);
l_z_um = vert_tiles(end).sample_xyz_um_done(3) + vert_tiles(end).stack_size(3) - ...
    ref_z_um;
merged_ll = [vert_tiles(1).stack_size(1:2), l_z_um];
merged_ll_um = [vert_tiles(1).stack_size_um(1:2), l_z_um];
ch_list = vert_tiles(1).channel;
num_ch = numel(ch_list);

tile_data = arrayfun(@(x) zeros(merged_ll_um, 'uint16'), 1:2, 'UniformOutput', false);
tile_im_1um = cell(num_ch, num_vt);
for i = 1 : num_vt
    tmp_im = vert_tiles(i).load_tile();
    tmp_z_mmxx = [vert_tiles(i).sample_xyz_um_done(3), vert_tiles(i).sample_xyz_um_done(3) + ...
        vert_tiles(i).stack_size_um(3) - 1] - ref_z_um + 1;
    for i_ch = ch_list
        tmp_im_rz = imresize3(medfilt3(tmp_im{i_ch}), vert_tiles(i).stack_size_um);
        tile_im_1um{i_ch, i} = tmp_im_rz;
        % Direct stitching
        tile_data{i_ch}(:, :, tmp_z_mmxx(1) : tmp_z_mmxx(2)) = ...
            max(tile_data{i_ch}(:, :, tmp_z_mmxx(1) : tmp_z_mmxx(2)), ...
            tmp_im_rz);
    end    
%     vert_tiles(i).clear_buffer();
end
%%
vis_im_yxz = fun_merge_image_stacks(tile_data);
vis_im_zyx = permute(vis_im_yxz, [4, 1, 3, 2]);
vis_im_zyx_mip = squeeze(max(vis_im_zyx(:, :, :, 150:250), [], 4));

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
image(ax_hdl, vis_im_zyx_mip);
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XLabel.String = 'Y (\mum)';
ax_hdl.YLabel.String = 'Z (\mum)';
%% Single xy tile z stitching - translation only 
vis_Q = false;
disp_vec = nan(3, num_vt, num_ch);
for reg_ch = 1 : num_ch
    for fix_idx = 1 : (num_vt - 1)
        im_fixed = tile_im_1um{reg_ch, fix_idx};
        im_move = tile_im_1um{reg_ch, fix_idx + 1};
        im_move(:, :, 1:20) = 0;
        pos_fixed_um = vert_tiles(fix_idx).sample_xyz_um_done;
        pos_mov_um = vert_tiles(fix_idx + 1).sample_xyz_um_done;
        pxl_shift = pos_mov_um - pos_fixed_um;
        pxl_shift(3) = pxl_shift(3);
        [pix_shift_e, max_c] = fun_im_reg_stack_masked_fft(im_fixed, im_move, 'pixel_shift_0', pxl_shift, ...
            'iadj', 3, 'masked_int', 1e4, 'vis_Q', vis_Q, 'search_range_yxz', [15, 15, 15]);
        disp_vec(:, fix_idx, reg_ch) = pix_shift_e;
        fprintf('Finish computing the displacement vector for tile %d channel %d\n', ...
            fix_idx, reg_ch);
    end
end
selected_disp_vec = cat(2, [0;0;0], disp_vec(:, 1:2, 2), disp_vec(:, 3:end-1, 1));
%% Construct image based on the estimated displacement vecotr 
stack_size_um = vert_tiles(1).stack_size_um;
cum_disp_vec = cumsum(selected_disp_vec, 2).';

disp_vec_sub_min = min(cum_disp_vec, [], 1);
disp_vec_sub_max = max(cum_disp_vec, [], 1);
ori_disp_vec = [0, 0, 0] - disp_vec_sub_min;
extended_im_size = [stack_size_um + disp_vec_sub_max - disp_vec_sub_min];
tile_data_fftr = arrayfun(@(x) zeros(extended_im_size, 'uint16'), ch_list, 'UniformOutput', false);
for ch = 1 : num_ch
    tmp_im_s = tile_data_fftr{ch};
    for i = 1 : num_vt
        tmp_im = tile_im_1um{ch, i};
%         tmp_im(:, :, 1:10) = 0;
        tmp_bbox_mmxx = cum_disp_vec(i, :) + ori_disp_vec;
        tmp_bbox_mmxx = [tmp_bbox_mmxx + 1, tmp_bbox_mmxx + stack_size_um];
        tmp_im_s(tmp_bbox_mmxx(1):tmp_bbox_mmxx(4), tmp_bbox_mmxx(2):tmp_bbox_mmxx(5), tmp_bbox_mmxx(3):tmp_bbox_mmxx(6)) = ...
            max(tmp_im_s(tmp_bbox_mmxx(1):tmp_bbox_mmxx(4), tmp_bbox_mmxx(2):tmp_bbox_mmxx(5), tmp_bbox_mmxx(3):tmp_bbox_mmxx(6)), ...
            tmp_im);
    end
    tile_data_fftr{ch} = tmp_im_s;
end

vis_imfr_yxz = fun_merge_image_stacks(tile_data_fftr);
vis_imfr_zyx = permute(vis_imfr_yxz, [4, 1, 3, 2]);
vis_imfr_zyx_mip = squeeze(max(vis_imfr_zyx(:, :, :, 150:250), [], 4));

avi_fp = fullfile(vis_folder, sprintf('%s_%s_single_tile_vstitch_grid_%d_zy.avi', ...
    exp_group, exp_name, vis_tile.grid_ind));
fun_vis_write_stack_to_avi(vis_imfr_zyx, avi_fp);
%%
fig_hdl = figure;
ax1 = subplot(1,2,1);
image(ax1, vis_im_zyx_mip);
ax1.DataAspectRatio = [1,1,1];
ax1.XLabel.String = 'Y (\mum)';
ax1.YLabel.String = 'Z (\mum)';
ax1.Title.String = 'Direct stitching';
ax2 = subplot(1,2,2);
image(ax2, vis_imfr_zyx_mip);
ax2.DataAspectRatio = [1,1,1];
ax2.XLabel.String = 'Y (\mum)';
ax2.YLabel.String = 'Z (\mum)';
ax2.Title.String = 'MFFT';
fig_fp = fullfile(vis_folder, sprintf('%s_%s_single_tile_vstitch_grid_%d_100_um_mip_zy.png', ...
    exp_group, exp_name, vis_tile.grid_ind));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
