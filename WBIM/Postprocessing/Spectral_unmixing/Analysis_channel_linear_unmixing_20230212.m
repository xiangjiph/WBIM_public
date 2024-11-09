DataManager = WBIMFileManager;
exp_group = 'TestSample';
exp_name = 'WBIM20220616002_20221027';
tile_str = DataManager.load_tile_in_experiment(exp_group, exp_name);
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'spectrum_unmixing');
%%
test_layer = 3;
test_tiles = tile_str.Scan{test_layer};
test_tiles = WBIMTileManager.select_duplicated_tiles(test_tiles);
[stitched_mip, local_tile_info] = WBIMTileManager.get_stitched_step_mip_1_ch(...
    test_tiles, 2);
[stitched_mip_1, local_tile_info_1] = WBIMTileManager.get_stitched_step_mip_1_ch(...
    test_tiles, 1);
%%
tv = WBIMTileViewer(test_tiles);
tv.update_tile_grid();
tv.init_step_mip(stitched_mip_1, local_tile_info_1.step_mip_pixel_yxz_um);
tv.ui_init_ctrl();
%% Load tiles and compute masks
num_layer = numel(tile_str.Scan);
layer_data = cell(num_layer, 1);
for i = 1 : num_layer
    tmp_tile = tile_str.Scan{i};
    tmp_data = struct;
    tmp_data.smip = cell(1, 2);
    for ch = 1 : 2
        [tmp_data.smip{ch}, tmp_data.local_tile_info] = WBIMTileManager.get_stitched_step_mip_1_ch(...
            tmp_tile, ch);
    end
    tmp_data.abl_str = WBIMAcqPostProcessing.compute_ablation_mask_yx_um_itp_vNs(...
        tmp_data.smip{1}, tmp_data.smip{2}, tmp_data.local_tile_info, 'visQ', false, ...
        'crmax_shgQ', false);
    layer_data{i} = tmp_data;
    fprintf('Finish computing the ablation mask for layer %d\n', i);
end
%% Try spectrum unmixing again
tt = test_tiles(73);
smip_cell = tt.load_step_mip();
smip_vsl = smip_cell{1};
smip_shg = smip_cell{2};
test_section = 2;
[stat_str, fig_hdl] = WBIMAcqPostProcessing.compute_spectrum_unmixing_factor(...
    smip_vsl(:, :, test_section), smip_shg(:, :, test_section), [5000, 32000], [0, 8000], true);
%%
umslope = nan(10, 5);
umr2 = nan(10, 5);
for test_layer = 1 : 5
    for test_section = 1 : 10
        im_vsl = layer_data{test_layer}.smip{1}(:, :, test_section);
        im_shg = layer_data{test_layer}.smip{2}(:, :, test_section);
        %%
        [stat_str, fig_hdl] = WBIMAcqPostProcessing.compute_spectrum_unmixing_factor(...
            im_vsl, im_shg, [1500, 32000], [0, 8000], false);
%         fig_fp = fullfile(vis_folder, 'spectrum_unmixing', sprintf('%s_%s_layer_%d_section_%d_spex_umix.png', ...
%             exp_group, exp_name, test_layer, test_section));
%         fun_print_image_in_several_formats(fig_hdl, fig_fp);
%         delete(fig_hdl);
        umslope(test_section, test_layer) = stat_str.k;
        umr2(test_section, test_layer) = stat_str.R2;
    end
end
%%
umslope = umslope(2:end, 2:end);
umr2 = umr2(2:end, 2:end);
ums_avg_r2 = mean(umslope .* umr2, 'all');
ums_avg2_r2 = mean(umslope .^2 .* umr2, 'all');
ums_std_r2 = sqrt(ums_avg2_r2 - ums_avg_r2.^2);
fig_hdl = figure;
ax_1 = nexttile();
scatter(ax_1, umslope(:), umr2(:));
ax_1.YLabel.String = 'R^2';
ax_1.XLabel.String = 'Slope';
legend(ax_1, sprintf('%.3f \\pm %.3f', ums_avg_r2, ums_std_r2), ...
    'Location', 'best');
ax_2 = nexttile();
imagesc(ax_2, umslope);
colorbar(ax_2);
ax_2.CLim = [0, 0.2];
ax_2.XLabel.String = 'Layer Index';
ax_2.YLabel.String = 'Section Index';
fig_fp = fullfile(vis_folder, 'spectrum_unmixing', sprintf('%s_%s_spex_umix_fitting_summary.png', ...
    exp_group, exp_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
% 