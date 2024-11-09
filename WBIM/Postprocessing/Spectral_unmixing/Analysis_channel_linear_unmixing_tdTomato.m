DataManager = WBIMFileManager;
exp_group = 'TestSample';
exp_name = 'WBIM20230414002_20230427';
% tile_str = DataManager.load_tile_in_experiment(exp_group, exp_name);
% vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
%     'spectrum_unmixing');
%%
im_folder = 'H:\WBIM\Acquisition\TestSample\WBIM20230414002_20230427\00006\Explore';
test_im_fp = fullfile(im_folder, 'WBIM_00005.tif');
test_im = DataManager.load_single_tiff(test_im_fp);
test_im_vsl = medfilt3(test_im(:, :, 1:2:end));
test_im_td = medfilt3(test_im(:, :, 2:2:end));
%%
vis_sec = 150;
fig_hdl = figure;
a1 = nexttile;
imagesc(a1, test_im_vsl(:, :, vis_sec));
a1.DataAspectRatio = [1,1,1];
a2 = nexttile;
imagesc(a2, test_im_td(:, :, vis_sec));
a2.DataAspectRatio = [1,1,1];
a3 = nexttile;
imagesc(a3,  max(0, test_im_td(:, :, vis_sec) - test_im_vsl(:, :, vis_sec)));
a3.DataAspectRatio = [1,1,1];
% histogram2(a3, test_im_vsl(:, :, vis_sec), test_im_td(:, :, vis_sec), ...
%     'DisplayStyle','tile');
a4 = nexttile;
imshowpair(test_im_vsl(:, :, vis_sec), test_im_td(:, :, vis_sec));




%% Try spectrum unmixing again
test_section = 150;
imshowpair(test_im_vsl(:, :, test_section), test_im_td(:, :, test_section));
[stat_str, fig_hdl] = WBIMAcqPostProcessing.compute_spectrum_unmixing_factor(...
    test_im_vsl(:, :, test_section), test_im_td(:, :, test_section),...
    [5000, 32000], [0, 8000], true);
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
%%

%% 3 channels spectral unmixing
% Implemented using WBIM20230126001_20230607
% Layer 1, exploration mode, tile list idx 36
sm_coeff = [  1, -0.12,      0, -0.12; ...
          -0.02,     1,      0, -0.02; ...
              0,     0,      1,     0; ...
              0,     0,      0,     1];
% Pick tiles from 
test_tile = process_tile(36);
test_tile.load_tile();
test_stacks = test_tile.tile;
test_stacks = cellfun(@(x) WBIMAcqPostProcessing.filter_sections(x, @medfilt2), test_stacks, 'UniformOutput', false);
% test_stacks = cellfun(@double, test_stacks, 'UniformOutput', false);

test_stacks_dm = WBIMAVD.smip_spectral_unmixing(test_stacks, 'coeff', sm_coeff);
    %% Visualization
vis_sec = 20;
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 3;
a1 = nexttile;
imagesc(a1, test_stacks{1}(:, :, vis_sec));
a1.Title.String = 'FITC';
a2 = nexttile;
imagesc(a2, test_stacks{2}(:, :, vis_sec));
a2.Title.String = 'SHG';
a3 = nexttile;
imagesc(a3, test_stacks{4}(:, :, vis_sec));
a3.Title.String = 'tdTomato';

a4 = nexttile; 
imagesc(a4, test_stacks_dm{1}(:, :, vis_sec));
a4.Title.String = 'FITC UM';
% imagesc(a4, max(0, test_stacks{3}(:, :, vis_sec) - 0.1 * test_stacks{1}(:, :, vis_sec)));
a5 = nexttile; 
imagesc(a5, test_stacks_dm{2}(:, :, vis_sec));
a5.Title.String = 'SHG UM';
% imagesc(a5, max(0, test_stacks{3}(:, :, vis_sec) - 0.1 * test_stacks{2}(:, :, vis_sec)));
a6 = nexttile; 
imagesc(a6, test_stacks_dm{4}(:, :, vis_sec));
a6.Title.String = 'tdTomato UM';
% imagesc(a6, max(0, test_stacks{3}(:, :, vis_sec) - 0.1 * test_stacks{1}(:, :, vis_sec) - 0.1 * test_stacks{2}(:, :, vis_sec)));

[a1.DataAspectRatio, a2.DataAspectRatio, a3.DataAspectRatio, ...
    a4.DataAspectRatio, a5.DataAspectRatio, a6.DataAspectRatio] = deal([1,1,1]);
a1.XLabel.String = 'X (\mum)';
a1.YLabel.String = 'Y (\mum)';
fig_fp = fullfile(vis_folder, sprintf('%s_%s_layer_%d_%s_tile_%d_sec_%d_unmixing.png', ...
    exp_group, exp_name, layer_idx, test_tile.acq_mode, test_tile.grid_ind, vis_sec));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
    %% Visualization 
figure;
imshowpair(test_stacks_dm{1}(:, :, vis_sec), test_stacks_dm{4}(:, :, vis_sec));
% imshowpair(test_stacks{4}(:, :, vis_sec), test_stacks{1}(:, :, vis_sec));