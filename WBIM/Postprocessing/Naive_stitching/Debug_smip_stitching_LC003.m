clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'Zhuhao';
exp_name = 'ZW2a3_20240518';
tiles_cell = DataManager.load_tile_in_experiment(exp_group, exp_name, WBIMMicroscopeMode.Scan);
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'smip_stitching');
%%
[smip_im, layer_data] = fun_visualization_stitch_smip_3D(tiles_cell(1:end-1), 'zero_last_frame_Q', true);
% Remove the last two frames: 
% smip_im = cellfun(@(x) x(:, :, 1:end-3), smip_im, 'UniformOutput', false);
% Flip the channel order: 
% smip_im{3} = smip_im{1};
% smip_im(1) = [];
%% Try spectral demising
% sd_smip = cell(size(smip_im));
% sd_smip{3} = rescale(single(smip_im{3} - smip_im{1} * 0.15) .^ 1/3);
%%
% Stretch contrast and merge channels
% rgb_im = fun_merge_image_stacks(smip_im, 'method', 'GraySkull', 'stretch_contrast_Q', true);
tile_data_sc = cell(numel(smip_im), 1);
tile_data_sc{1} = im2uint8(rescale(single(smip_im{1})).^(1/2));
tile_data_sc{2} = im2uint8(rescale(single(smip_im{2})).^(1/2));
tile_data_sc{3} = im2uint8(rescale(single(smip_im{3})).^(1/2));
% tmp = tile_data_sc{1};
% tile_data_sc{1} = tile_data_sc{2};
% tile_data_sc{2} = tmp;

rgb_im = fun_merge_image_stacks(tile_data_sc, 'method', 'GraySkull', 'stretch_contrast_Q', false);
rgb_im_zyx = permute(rgb_im, [4,1,3,2]);
rgb_im_zyx = rgb_im_zyx(:, :, :, end:-1:1); % start from the front. 
rgb_im_zyx = rgb_im_zyx(:, end:-1:1, :, :); % flip to match Zhuhao's visualization
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
layer_idx = 1;
[smip_im, layer_data] = fun_visualization_stitch_smip_3D(tiles_cell(layer_idx), 'zero_last_frame_Q', true);
tv = WBIMTileViewer(tiles_cell{layer_idx});
tv.update_tile_grid();
tv.init_step_mip(layer_data(layer_idx).smip{3}, layer_data(layer_idx).local_tile_info.step_mip_pixel_yxz_um);
tv.ui_init_ctrl();
% %%
% figure;imagesc(shg_mip);set(gca, 'DataAspectRatio', [1,1,1]);
%% 
num_layer = numel(layer_data);
layer_smip = cell(num_layer, 1);
for i = 1 : num_layer
    tmp_data = layer_data(i).smip;
    tmp_data = cellfun(@(x) im2uint8(rescale(x).^(1/2)), tmp_data, 'UniformOutput', false);
    layer_smip{i} = fun_merge_image_stacks(tmp_data, 'method', 'GraySkull', 'stretch_contrast_Q', false);
end
%%
imfilter