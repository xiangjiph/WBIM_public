clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'TestSample';
exp_name_normal = 'WBIM20230125001_20230213';
exp_name_decal = 'WBIM20230124001_20230301';
tile_normal = DataManager.load_tile_in_experiment(exp_group, exp_name_normal);
tile_decal = DataManager.load_tile_in_experiment(exp_group, exp_name_decal);
exp_tiles = {tile_normal.Scan, tile_decal.Scan};
%%
layer_idx = 1;

tile_cell = {tile_normal.Scan{layer_idx}, tile_decal.Scan{layer_idx}};
    
[smip_s_n, region_info_n] = WBIMTileManager.get_stitched_step_mip_1_ch(...
    tile_cell{1}, 2);

[smip_s_d, region_info_d] = WBIMTileManager.get_stitched_step_mip_1_ch(...
    tile_cell{2}, 2);

%%
tv_n = WBIMTileViewer(tile_cell{1});
tv_n.update_tile_grid();
tv_n.init_step_mip(smip_s_n,region_info_n.step_mip_pixel_yxz_um);
% tv.init_step_mip_mask(new_abl_str.mask{2});
% tv.init_step_mip_mask(vsl_result.abl_mask);
tv_n.ui_init_ctrl();
%%
tv_d = WBIMTileViewer(tile_cell{2});
tv_d.update_tile_grid();
tv_d.init_step_mip(smip_s_d,region_info_d.step_mip_pixel_yxz_um);
tv_d.ui_init_ctrl();
%%
test_exp = 2;
test_tile = tile_cell{test_exp}(8);

test_stacks = test_tile.load_tile([1, 2]);
test_stacks = cellfun(@medfilt3, test_stacks, 'UniformOutput', false);
% Downsample the image stack
test_stacks = cellfun(@(x) imresize3(x, round(test_tile.stack_size_um)),...
    test_stacks, 'UniformOutput', false);
test_mips = cellfun(@(x) max(x, [], 3), test_stacks, 'UniformOutput', false);
skull_pixels = (imerode(test_mips{1} < 6e3, strel('disk', 10))) & (test_mips{2} > 5e3);
skull_ind = find(skull_pixels);
num_skull_ind = numel(skull_ind);
shg_z_profile = reshape(permute(test_stacks{2}, [3,1,2]), test_tile.stack_size(3), []);
shg_z_profile = shg_z_profile(:, skull_ind);


[ls_um, ls_fit_R2] = deal(nan(num_skull_ind, 1));
profile_z_um = (1 : test_tile.stack_size(3));
est_bg_max = 3e3;
t_tic = tic;
parfor i = 1 : num_skull_ind
    warning('off', 'stats:LinearModel:RankDefDesignMat');
    tmp_int = shg_z_profile(:, i);
    if any(est_bg_max < tmp_int)
        tmp_str = WBIMAcqPostProcessing.compute_z_int_scatter_length_from_exp_fit(...
            tmp_int, profile_z_um, 'visQ', false, 'background_int', est_bg_max, ...
            'saturation_int', 3.1e4);
        ls_um(i) = tmp_str.ls_um;
        ls_fit_R2(i) = tmp_str.R2;
    end
    if mod(i, 1000) == 0
       fprintf('Finish computing %d/%d profiles. Elapsed time is %.2f seconds\n',...
           i, num_skull_ind, toc(t_tic)); 
    end
    warning('on', 'stats:LinearModel:RankDefDesignMat');
end
%%
ls_map = nan(size(test_mips{1}));
is_valid_ls_Q = (ls_um > 10) & (abs(ls_um) < 1e2);
ls_map(skull_ind(is_valid_ls_Q)) = ls_um(is_valid_ls_Q);
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
ax_1 = subplot(1,2,1);
histogram2(ax_1, ls_um(is_valid_ls_Q), ls_fit_R2(is_valid_ls_Q), 'DisplayStyle', 'tile');
med_ls_um = median(ls_um(is_valid_ls_Q));
ax_1.Title.String = sprintf('Median: %.1f \\mum', med_ls_um);
ax_1.XLabel.String = 'Scattering length (\mum)';
ax_1.YLabel.String = 'R^2';
cbar = colorbar(ax_1);
cbar.Label.String = 'Number of points';
ax_1.YLim(1) = 0.80;
ax_1.ColorScale = 'linear';
ax_2 = subplot(1,2,2);
imagesc(ax_2, ls_map);
ax_2.DataAspectRatio = [1,1,1];
cbar_hdl = colorbar(ax_2);
ax_2.CLim = [0, 100];
cbar_hdl.Label.String = 'l_s (\mum)';

fig_fp = fullfile(DataManager.fp_experiment(exp_group, test_tile.experiment),...
    'visualization', 'Scattering_length', ...
    sprintf('%s_%s_layer_%d_grid_%d_scatter_length_map.png', ...
    exp_group, test_tile.experiment, test_tile.layer, test_tile.grid_ind));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Processing layer
test_exp = 2;
test_layer = 2;
layer_tiles = exp_tiles{test_exp}{test_layer};
[smip_s_d, region_info_d] = WBIMTileManager.get_stitched_step_mip_1_ch(...
    layer_tiles, 2);
tv_d = WBIMTileViewer(layer_tiles);
tv_d.update_tile_grid();
tv_d.init_step_mip(smip_s_d,region_info_d.step_mip_pixel_yxz_um);
tv_d.ui_init_ctrl();
%%
tile_ls_stat = WBIMAcqPostProcessing.analyze_stack_scattering_length_in_skull(layer_tiles(40), ...
    'visQ', true, 'save_fig_Q', true);