clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'TestSample';
exp_name = 'WBIM20230425001_20230515';
% DataManager.download_tile_info_in_experiment(exp_group, exp_name);
tile_str = DataManager.load_tile_in_experiment(exp_group, exp_name, WBIMMicroscopeMode.Scan);
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'smip_stitching');
%% Load tiles SMIP
process_tiles = tile_str;
valid_layer_Q = ~cellfun(@isempty, process_tiles);
process_tiles = process_tiles(valid_layer_Q);
process_tiles_s = cellfun(@(x) WBIMTileManager.select_duplicated_tiles(x), process_tiles, 'UniformOutput', false);
% num_layer = numel(process_tiles_s);
num_layer = 4;
ch_list = process_tiles_s{1}(1).channel;
num_ch = numel(ch_list);
layer_data = cell(num_layer, 1);
for i = 1 : num_layer
    tmp_tile = process_tiles_s{i};
    tmp_data = struct;
    tmp_data.smip = cell(1, num_ch);
    for ch = 1 : num_ch
        [tmp_data.smip{ch}, tmp_data.local_tile_info] = WBIMTileManager.get_stitched_step_mip_1_ch(...
            tmp_tile, ch_list(ch));
    end
    arrayfun(@(x) x.clear_buffer(), tmp_tile);
    layer_data{i} = tmp_data;
    fprintf('Finish loading data in layer %d\n', i);
end
layer_data = cat(1, layer_data{:});
%%
smip_pixel_yxz_um = layer_data(1).local_tile_info.step_mip_pixel_yxz_um;
stack_size = layer_data(1).local_tile_info.stack_size;
vol_bbox_xy_mm_um = arrayfun(@(x) x.local_tile_info.region_bbox_mm_um, layer_data, 'UniformOutput', false);
vol_bbox_xy_mm_um = cat(1, vol_bbox_xy_mm_um{:});
vol_bbox_z_mm_um = arrayfun(@(x) x.local_tile_info.stage_z_um_done , layer_data);
vol_bbox_xy_xx_um = arrayfun(@(x) x.local_tile_info.region_bbox_xx_um, layer_data, 'UniformOutput', false);
vol_bbox_xy_xx_um = cat(1, vol_bbox_xy_xx_um{:});

vol_bbox_mmxx_um = cat(2, vol_bbox_xy_mm_um, vol_bbox_z_mm_um, vol_bbox_xy_xx_um, vol_bbox_z_mm_um + stack_size(3) - 1);
vol_bbox_mmxx_sp = round(vol_bbox_mmxx_um ./ [smip_pixel_yxz_um, smip_pixel_yxz_um]);
space_mm_sp = min(vol_bbox_mmxx_sp(:, 1:3), [], 1);
space_xx_sp = max(vol_bbox_mmxx_sp(:, 4:6), [], 1);
space_ll_sp = space_xx_sp - space_mm_sp + 1;
vol_bbox_mmxx_sp_l = vol_bbox_mmxx_sp - [space_mm_sp, space_mm_sp] + 1;
%%
stitch_im_type = 'uint16';
num_s_ch = numel(layer_data(1).smip);
smip_im = cell(num_s_ch, 1);
for i_c = 1 : num_s_ch
    tmp_smip_im = zeros(space_ll_sp, stitch_im_type);
    for i_l = 1 : num_layer
        tmp_im = layer_data(i_l).smip{i_c};
        tmp_im_size = size(tmp_im);
        tmp_bbox_mm = vol_bbox_mmxx_sp_l(i_l, 1:3);
        tmp_bbox_xx = tmp_bbox_mm + tmp_im_size - 1;
        tmp_smip_im(tmp_bbox_mm(1) : tmp_bbox_xx(1), tmp_bbox_mm(2) : tmp_bbox_xx(2), tmp_bbox_mm(3) : tmp_bbox_xx(3)) = ...
            max(tmp_smip_im(tmp_bbox_mm(1) : tmp_bbox_xx(1), tmp_bbox_mm(2) : tmp_bbox_xx(2), tmp_bbox_mm(3) : tmp_bbox_xx(3)), uint16(tmp_im));        
        fprintf('Finish adding layer %d of channel %d\n', i_l, i_c);
    end        
    smip_im{i_c} = tmp_smip_im;
end
%% Try spectral demising
sd_smip = cell(size(smip_im));
sd_smip{3} = rescale(single(smip_im{3} - smip_im{1} * 0.15) .^ 1/3);

%%
% Stretch contrast and merge channels
% rgb_im = fun_merge_image_stacks(smip_im, 'method', 'GraySkull', 'stretch_contrast_Q', true);
tile_data_sc = cell(num_ch, 1);
tile_data_sc{1} = im2uint8(rescale(single(smip_im{1})).^(1/2));
tile_data_sc{2} = im2uint8(rescale(single(smip_im{2})).^(1/3));
% tile_data_sc{3} = im2uint8(rescale(single(smip_im{3})).^(1/2));
rgb_im = fun_merge_image_stacks(tile_data_sc, 'method', 'GraySkull', 'stretch_contrast_Q', false);
rgb_im_zyx = permute(rgb_im, [4,1,3,2]);
rgb_im_zxy = permute(rgb_im, [4,2,3,1]);

avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_merged_ds5_yxz_sc.avi', ...
    exp_group, exp_name));
fun_vis_write_stack_to_avi(rgb_im, avi_fp, 'Quality', 90);
cvrt_out = fun_vis_avi_to_mp4(avi_fp);

DataManager.write_tiff_stack(rgb_im, strrep(avi_fp, 'avi', 'tif'), 'color');

avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_merged_ds5_zyx_sc.avi', ...
    exp_group, exp_name));
fun_vis_write_stack_to_avi(rgb_im_zyx, avi_fp, 'Quality', 90);
cvrt_out = fun_vis_avi_to_mp4(avi_fp);

avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_merged_ds5_zxy_sc.avi', ...
    exp_group, exp_name));
fun_vis_write_stack_to_avi(rgb_im_zxy, avi_fp, 'Quality', 90);
cvrt_out = fun_vis_avi_to_mp4(avi_fp);
%%
% tv = WBIMTileViewer(tile_str.Scan{layer_idx});
% tv.update_tile_grid();
% tv.init_step_mip(layer_data{layer_idx}.smip{2}, layer_data{layer_idx}.local_tile_info.step_mip_pixel_yxz_um);
% tv.ui_init_ctrl();
% %%
% figure;imagesc(shg_mip);set(gca, 'DataAspectRatio', [1,1,1]);













