% This script is for generating the dataset for testing Janelia's rendering
% software. 
clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'Zhuhao';
exp_name = 'ZW2a1_20240502';

tiles_cell = DataManager.load_tile_in_experiment(exp_group, exp_name, WBIMMicroscopeMode.Scan);

all_scan_tile = cellfun(@(x) WBIMTileManager.select_duplicated_tiles(x, 'latest'), ...
    tiles_cell, 'UniformOutput', false);
all_scan_tile = cat(1, all_scan_tile{:});
volume_info = WBIMTileMetadata.combine_tile_info_3D(all_scan_tile);
%% Load explored and scanned tiles
layer_idx = 4;
ch_list = [1];
process_tile = tiles_cell{layer_idx};
process_tile = WBIMTileManager.select_duplicated_tiles(process_tile,'latest');
%%
[smip_cells, local_tile_info] = WBIMTileManager.get_stitched_step_mip(process_tile, ...
    ch_list);
%% single layer visualization
tv = WBIMTileViewer(process_tile);
tv.update_tile_grid();
% tv.init_step_mip(smip_cells{4}, local_tile_info.step_mip_pixel_yxz_um);
tv.init_step_mip(smip_cells{1}, local_tile_info.step_mip_pixel_yxz_um);
% tv.init_step_mip_mask(avd_hdl.mask_cell{2});
% tv.init_step_mip_mask(vsl_result.abl_mask);
tv.ui_init_ctrl();
%% Figure out when the backgroud increase dramatically
% Explore layer 9: background level increase dramatically at tile 136
layer_acq_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), ...
    '00009_disruppted', char(WBIMMicroscopeMode.Explore));
tile_dir = dir(fullfile(layer_acq_folder, '**', DataManager.fn_acq_info));
num_tiles = numel(tile_dir);
tile_list = cell(num_tiles, 1);
if num_tiles > 0
    for i = 1 : num_tiles
        tile_list{i} = WBIMTileMetadata.load(fullfile(tile_dir(i).folder,...
            tile_dir(i).name));
    end
    tile_list = cat(1, tile_list{:});
end
tile_done_time = arrayfun(@(x) x.t_done, tile_list, 'UniformOutput', false);
[s_t, s_ind] = sort(tile_done_time);
tile_list_1 = tile_list(s_ind).copy();
tile_dir = tile_dir(s_ind);
num_tiles = numel(tile_list_1);
check_ch = 4;
mip_fp_cell = arrayfun(@(x) fullfile(x.folder, DataManager.fn_acq_mip(check_ch)), ...
    tile_dir, 'UniformOutput', false);
tile_mip_cell = cellfun(@(x) DataManager.load_single_tiff(x), mip_fp_cell, ...
    'UniformOutput', false);
% tile_mip_cell = cat(1, tile_mip_cell{:});

local_im = Grid.stitch_tiles(cat(1, tile_list_1.tile_mmxx_pxl), tile_mip_cell, 'mip');   
stat = struct;
[stat.avg_int, stat.std_int, stat.max_adj_diff_int] = deal(nan(num_tiles, 1));
for i = 1 : num_tiles
    tmp_tile = single(tile_mip_cell{i});
    tmp_bg_int = tmp_tile(tmp_tile < 6000);
    stat.avg_int(i) = mean(tmp_bg_int);
    stat.std_int(i) = std(tmp_bg_int);    
end
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
plot(stat.avg_int);
hold(ax_hdl, 'on');
plot(stat.std_int);
%%
tile_idx = 135;
check_tile = tile_list_1(tile_idx);
check_stat = check_tile.load_stat();
check_mip = check_tile.load_mip();
check_mip_merge = fun_merge_image_stacks(check_mip([2, 1, 3]), 'method', 'GraySkull');
% check_smip = check_tile.load_step_mip();
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
image(ax_hdl, check_mip_merge);
ax_hdl.Title.String = num2str(tile_idx);
fprintf('******\n');
fprintf('Tile %d acquired at %s\n', tile_idx, check_tile.t_done);
disp(check_stat.CH4.z_stat.mean');
%%
vis_ch = 3;
[smip_cells, local_tile_info] = WBIMTileManager.get_stitched_step_mip(tile_list_1, ...
    vis_ch);
tile_mip_4 = WBIMTileManager.get_stitched_tile_mip(tile_list_1, 4);
tile_mip_3 = WBIMTileManager.get_stitched_tile_mip(tile_list_1, 3);
tile_mip_2 = WBIMTileManager.get_stitched_tile_mip(tile_list_1, 2);

tv = WBIMTileViewer(tile_list_1);
tv.update_tile_grid();
tv.init_step_mip(smip_cells{vis_ch}, local_tile_info.step_mip_pixel_yxz_um);
tv.ui_init_ctrl();
%% 05/08/2024 during scan
tile_list = DataManager.load_tile_in_layer(exp_group, exp_name, 12, WBIMMicroscopeMode.Scan);
tile_done_time = arrayfun(@(x) x.t_done, tile_list, 'UniformOutput', false);
[s_t, s_ind] = sort(tile_done_time);
tile_list = tile_list(s_ind);

vis_ch = 3;
[smip_cells, local_tile_info] = WBIMTileManager.get_stitched_step_mip(tile_list, ...
    vis_ch);
tv = WBIMTileViewer(tile_list);
tv.update_tile_grid();
tv.init_step_mip(smip_cells{vis_ch}, local_tile_info.step_mip_pixel_yxz_um);
tv.ui_init_ctrl();

