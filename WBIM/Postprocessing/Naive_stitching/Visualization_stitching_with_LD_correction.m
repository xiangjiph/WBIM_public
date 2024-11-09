% Stitching with static lens deformation correction 
clc;clear;close all;
DataManager = WBIMFileManager;
%%
exp_group = 'Ablation_test';
exp_name = 'WBIMSkl_20240305002_20240709_3';
tile_str = DataManager.load_tile_in_experiment(exp_group, exp_name);
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'Stitched');
%% Parameters
write_stitched_data_Q = false;
correct_ld_Q = false;

stitch_voxel_size_um = [1, 1, 1];
zero_num_sec = 0;
zero_last_section_Q = true;
medfilt_Q = false;
%% Test deformation code on a single tile
if correct_ld_Q
    def_exp_name = '20230403_fiber';
    def_folder_root = DataManager.fp_layer('Calibration', def_exp_name, 1);
    def_mode = char(WBIMMicroscopeMode.Scan);
    def_fp = fullfile(def_folder_root, def_mode, 'polyCoeffVars.mat');
    assert(isfile(def_fp), 'Deformation file does not exist');
    def_str = DataManager.load_data(def_fp);
    
    test_tile = tile_str.(def_mode){10}(40);
    %Calculate polynomial deformation points
    [~, ~, gridSamplePts] =  invertPolyTransform(def_str.kBeta, def_str.d, ...
        test_tile.stack_size, true);
    
    im_stacks = test_tile.load_tile();
    
    tic
    transI = deformImageTiles(im_stacks{1}, gridSamplePts, 'flyback_Q', false, ...
        'interpMethod', 'gridded', 'processMethod', 'volume');
    toc
end
%% Load all the tiles 
t_tic = tic;
% stitch_set = WBIMMicroscopeMode.Scan;
% stitch_tiles = tile_str.(char(stitch_set));
% layer_list = 1:4;
stitch_tiles = scan_tiles_cell(2:2:end);
layer_list = 1 : numel(stitch_tiles);
stitch_tiles = cat(1, stitch_tiles{layer_list});
ch_list = stitch_tiles(1).channel;
num_ch = numel(ch_list);
num_tiles = numel(stitch_tiles);
% Compute overall bounding box
bbox_mmxx_um = cat(1, stitch_tiles.tile_mmxx_um);
layer_z_um = [stitch_tiles.layer_z_um];
stack_size_um = stitch_tiles(1).stack_size_um;

stack_size = stitch_tiles(1).stack_size;
ds_stack_size = round(stack_size_um ./ stitch_voxel_size_um);

bbox_mmxx_pxl = cat(1, stitch_tiles.tile_mmxx_pxl);

vol_bbox_z_mx_um = [min(layer_z_um), max(layer_z_um) + stack_size_um(3) - 1];
vol_bbox_mm_um = [min(bbox_mmxx_um(:, 1:2), [], 1), vol_bbox_z_mx_um(1)];
vol_bbox_xx_um = [max(bbox_mmxx_um(:, 3:4), [], 1), vol_bbox_z_mx_um(2)];
vol_bbox_ll_um = vol_bbox_xx_um - vol_bbox_mm_um + 1;

ds_bbox_ll = round(vol_bbox_ll_um ./ stitch_voxel_size_um);

tile_data = cell(1, num_ch);
% This process is dominated by reading H5 file (0.3 second per tile)
for i_ch = 1 : num_ch
    tmp_ch = ch_list(i_ch);
    tmp_stitch_data = zeros(ds_bbox_ll, 'uint16');
    for i = 1 : num_tiles
        try
            tmp_tile = stitch_tiles(i);
            tmp_tile_data = tmp_tile.load_tile(tmp_ch);
            tmp_tile_data = tmp_tile_data{1};
            % Apply lens deformation correction 
            if correct_ld_Q
                tmp_tile_data = deformImageTiles(tmp_tile_data, gridSamplePts, 'flyback_Q', false, ...
                    'interpMethod', 'gridded', 'processMethod', 'volume');
            end
            if medfilt_Q
                tmp_tile_data = medfilt3(tmp_tile_data);
            end
            if zero_num_sec && (tmp_tile.layer > 1)
                tmp_tile_data(:, :, 1:zero_num_sec) = 0;
            end
            if zero_last_section_Q
                tmp_tile_data(:, :, end) = 0;
            end
            tmp_tile_bbox_mm_um = [tmp_tile.tile_mmxx_um(1:2), tmp_tile.layer_z_um];
            tmp_tile_bbox_ll_um = [tmp_tile.tile_mmll_um(3:4), tmp_tile.stack_size_um(3)];
            tmp_tile_ll_ds_pxl = round(tmp_tile_bbox_ll_um ./ stitch_voxel_size_um);
            % Downsample image stack - need smoothing? 
            tmp_tile_data = imresize3(tmp_tile_data, tmp_tile_ll_ds_pxl);
            tmp_tile.clear_buffer();
            % Local bounding box 
            tmp_local_bbox_um = tmp_tile_bbox_mm_um - vol_bbox_mm_um;
            tmp_local_bbox_mm_ds_pxl = round(tmp_local_bbox_um ./ stitch_voxel_size_um);
            % Deal with edge: 
            tmp_local_bbox_mm_ds_pxl = max(tmp_local_bbox_mm_ds_pxl, 1);
            tmp_local_bbox_xx_ds_pxl = tmp_local_bbox_mm_ds_pxl + tmp_tile_ll_ds_pxl - 1;
            % Max - rendering
            tmp_stitch_data(tmp_local_bbox_mm_ds_pxl(1) : tmp_local_bbox_xx_ds_pxl(1), ...
                tmp_local_bbox_mm_ds_pxl(2) : tmp_local_bbox_xx_ds_pxl(2), ...
                tmp_local_bbox_mm_ds_pxl(3) : tmp_local_bbox_xx_ds_pxl(3)) = max(tmp_stitch_data(...
                tmp_local_bbox_mm_ds_pxl(1) : tmp_local_bbox_xx_ds_pxl(1), ...
                tmp_local_bbox_mm_ds_pxl(2) : tmp_local_bbox_xx_ds_pxl(2), ...
                tmp_local_bbox_mm_ds_pxl(3) : tmp_local_bbox_xx_ds_pxl(3)), tmp_tile_data);
            fprintf('Finish adding tile %d (%.3f %%)\n', i, (i/num_tiles) * 100);
        catch ME
            fprintf('Failed to add tile %d (%.3f %%)\n', i, (i/num_tiles) * 100);
        end
    end
    tile_data{i_ch} = tmp_stitch_data;
    % Write to tiff stack? 
    if write_stitched_data_Q
        stack_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_CH_%d.tif', ...
            exp_group, exp_name, i_ch));
        DataManager.write_tiff_stack(tile_data{i_ch}, stack_fp);
    end
    fprintf('Finish processing channel %d. Elapsed time is %.2f seconds\n', ...
        i_ch, toc(t_tic));
end
%%
i_ch = 2;
avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_CH_%d.avi', ...
            exp_group, exp_name, i_ch));
fun_vis_write_stack_to_avi(rescale(single(tile_data{i_ch})), avi_fp);
%% Merge channel 
rgb_im = zeros([ds_bbox_ll, 3], 'uint8');
rgb_im(:, :, :, 1) = im2uint8(fun_stretch_contrast(tile_data{1}));
rgb_im(:, :, :, 2) = im2uint8(fun_stretch_contrast(tile_data{2}));
rgb_im = permute(rgb_im, [1,2,4,3]);
% avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_merged.avi', ...
%     exp_group, exp_name));
% fun_vis_write_stack_to_avi(rgb_im, avi_fp);
%%
rgb_im_zyx = permute(rgb_im, [3, 1, 4, 2]);
% implay(rgb_im_zyx);
avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_merged_zyx.avi', ...
    exp_group, exp_name));
fun_vis_write_stack_to_avi(rgb_im_zyx, avi_fp);

% rgb_im = fun_merge_image_stacks(tile_data, 'method', 'GraySkull', 'stretch_contrast_Q', true);
%% Downsample 2x 
% Convert to log
% rgb_im = fun_merge_image_stacks(tile_data, 'method', 'rgb', ...
%     'stretch_contrast_Q', true);

tile_data_sc = cell(num_ch, 1);
% Increase by 1 before taking the log 
% tile_data_sc{1} = im2uint8(fun_stretch_contrast(tile_data{1}));
tile_data_sc{1} = im2uint8(rescale(single(tile_data{1})).^(1/2));
tile_data_sc{2} = im2uint8(rescale(single(tile_data{2})).^(1/3));
tile_data_sc{3} = im2uint8(rescale(single(tile_data{3})).^(1/3));
% for i = 1 : num_ch
%     tile_data_sc{i} = im2uint8(rescale(single(tile_data{i})).^(1/2));
% end
rgb_im = fun_merge_image_stacks(tile_data_sc, 'method', 'GraySkull', 'stretch_contrast_Q', false);
% rgb_im = cat(4, tile_data_sc{[1,3,2]});
% rgb_im = permute(rgb_im, [1,2,4,3]);


% tile_data_sc = cellfun(@fun_stretch_contrast, tile_data, 'UniformOutput', false);
% tmp = fun_stretch_contrast(tile_data{3});
% tmp_zyx = permute(tmp, [3, 1, 2]);
% implay(tmp_zyz);
% ds_ds_bbox_ll = round(ds_bbox_ll / 2);
% im_ch1_ds2 = im2uint8(fun_stretch_contrast(imresize3(tile_data{1}, ds_ds_bbox_ll)));
% im_ch2_ds2 = im2uint8(fun_stretch_contrast(imresize3(tile_data{2}, ds_ds_bbox_ll)));
% im_ch3_ds2 = im2uint8(fun_stretch_contrast(imresize3(tile_data{3}, ds_ds_bbox_ll)));
% rbg_im = repmat(im_ch2_ds2, 1, 1, 1, 3) * 0.6;
% rbg_im(:, :, :, 1) = max(im_ch1_ds2, rbg_im(:, :, :, 1));
% rbg_im(:, :, :, 2) = max(im_ch3_ds2, rbg_im(:, :, :, 2));
% 
% rbg_im = permute(rbg_im, [1,2,4,3]);

rgb_im_zyx = permute(rgb_im, [4,1,3,2]);
rgb_im_zxy = permute(rgb_im, [4,2,3,1]);

avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_merged_ds2_yxz_sc.avi', ...
    exp_group, exp_name));
fun_vis_write_stack_to_avi(rgb_im, avi_fp);
% Conver to MP4

avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_merged_ds2_zyx_sc.avi', ...
    exp_group, exp_name));
fun_vis_write_stack_to_avi(rgb_im_zyx, avi_fp);
DataManager.write_tiff_stack(rgb_im_zyx, strrep(avi_fp, 'avi', 'tif'), 'color');

avi_fp = fullfile(vis_folder, sprintf('%s_%s_stitched_stack_merged_ds2_zxy_sc.avi', ...
    exp_group, exp_name));
fun_vis_write_stack_to_avi(rgb_im_zxy, avi_fp);
DataManager.write_tiff_stack(rgb_im_zxy, strrep(avi_fp, 'avi', 'tif'), 'color');
%%
vis_mip = squeeze(max(rgb_im(:, :, :, 1:100), [], 4));

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 4;
ax_hdl = axes(fig_hdl);
image(ax_hdl, vis_mip);
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.Title.String = sprintf('Top %d \\mum MIP', 100 * 2);
ax_hdl.XLabel.String = 'X/2\mum';
ax_hdl.YLabel.String = 'Y/2\mum';
fig_fp = fullfile(vis_folder, sprintf('%s_%s_top_200_um_mip.png', exp_group,...
    exp_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
start_idx = 1400;
num_sec = 50;
vis_mip_zyx = squeeze(max(rgb_im_zyx(:, :, :, start_idx:(start_idx + num_sec)), [], 4));
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 4;
ax_hdl = axes(fig_hdl);
image(ax_hdl, vis_mip_zyx);
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.Title.String = sprintf('%d \\mum MIP', 50 * 2);
ax_hdl.XLabel.String = 'Y/2\mum';
ax_hdl.YLabel.String = 'Z/2\mum';
fig_fp = fullfile(vis_folder, sprintf('%s_%s_zyx_%d_%d_mip.png', exp_group,...
    exp_name, start_idx, start_idx + num_sec));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

%%
mip_range = 1700 : 1950;
mip_ch1 = max(tile_data{1}(:, :, mip_range), [], 3) * 2;
mip_ch2 = max(tile_data{2}(:, :, mip_range), [], 3) * 2;
mip_rgb = repmat(mip_ch2, 1, 1, 3) * 0.6;
mip_rgb(:, :, 1) = max(mip_rgb(:, :, 1), mip_ch1);

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
imshow(mip_rgb);
ax_hdl.Visible = 'on';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
fig_fp = fullfile(vis_folder, sprintf('%s_%s_vsl_mip_%d_to_%d_um.png', exp_group, exp_name, ...
    mip_range(1), mip_range(end)));
fun_print_image(fig_hdl, fig_fp);
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 3;
ax_hdl = axes(fig_hdl);
im_hdl = image(ax_hdl, mip_rgb);
ax_hdl.DataAspectRatio = [1,1,1];
% colormap(ax_hdl, 'gray');
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
% colormap('gray');
fig_fp = fullfile(vis_folder, sprintf('%s_%s_vsl_mip_%d_to_%d_um.png', exp_group, exp_name, ...
    mip_range(1), mip_range(2)));
fun_print_image(fig_hdl, fig_fp);
%% ZY view
% mip_range = 500 : 550;
mip_range = 560 : 595;
zy_mip = max(rg_im_zxy(:, :, :, mip_range), [], 4);

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
image(ax_hdl, zy_mip);
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XLabel.String = 'Right - Left (\mum)';
ax_hdl.YLabel.String = 'Ventral - Dorsal (\mum)';

ax_hdl = fun_set_axis_ticklabel_with_pixel_size(ax_hdl, 2);
ax_hdl.Title.String = sprintf('MIP %d \\mum to %d \\mum', mip_range(1) * 2 , 2 * mip_range(end));
fig_fp = fullfile(vis_folder, sprintf('%s_%s_vsl_mip_zy_%d_to_%d_um.png', exp_group, exp_name, ...
    mip_range(1) * 2, mip_range(end) * 2));
fun_print_image(fig_hdl, fig_fp);
%% ZX view
mip_range = 1450: 1500;
zx_mip = max(rgb_im_zxy(:, :, :, mip_range), [], 4);

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
image(ax_hdl, zx_mip );
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XTickLabel = arrayfun(@(x)num2str(2 * x), ax_hdl.XTick, 'UniformOutput', false);
ax_hdl.YLabel.String = 'Ventral - Dorsal (\mum)';
ax_hdl.XLabel.String = 'Posterior - Anterior (\mum)';
ax_hdl = fun_set_axis_ticklabel_with_pixel_size(ax_hdl, 2);
ax_hdl.Title.String = sprintf('MIP %d \\mum to %d \\mum', mip_range(1) * 2 , 2 * mip_range(end));
fig_fp = fullfile(vis_folder, sprintf('%s_%s_vsl_mip_zx_%d_to_%d_um.jpg', exp_group, exp_name, ...
    mip_range(1) * 2, mip_range(end) * 2));
fun_print_image(fig_hdl, fig_fp);
