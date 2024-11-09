vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'Analysis');
%% Compute the acquisition duration 
acq_time_fun = @(x) seconds(datetime(x.t_done) - datetime(x.t_init));
acq_time = cell(size(tiles_cell));
num_layer = numel(acq_time);
for i = 1 : numel(acq_time)
    tmp_layer_tiles = tiles_cell{i};
    acq_time{i} = arrayfun(acq_time_fun, tmp_layer_tiles);
end
% Remove empty layer
layer_idx = 1 : num_layer;
is_valid_Q = ~cellfun(@isempty, acq_time);
layer_idx = layer_idx(is_valid_Q);
acq_time = acq_time(is_valid_Q);
% Get statistics
acq_time_stat = cellfun(@fun_analysis_get_basic_statistics, acq_time, 'UniformOutput', false);
%%
acq_time_stat = cat(1, acq_time_stat{:});
acq_time_prtl = arrayfun(@(x) x.prctile_val([7, 8, 9]), acq_time_stat, 'UniformOutput', false);
acq_time_prtl = cat(1, acq_time_prtl{:});
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
p_m = plot(ax_hdl, layer_idx, acq_time_prtl(:, 2), 'LineWidth', 1.5);
p_c = fun_vis_confidence_interval_shaded(layer_idx, acq_time_prtl(:, 2), ...
    acq_time_prtl(:, 1), acq_time_prtl(:, 3), ax_hdl);
%% Box plot
bp_y = cat(1, acq_time{:});
bp_x = repelem(layer_idx, cellfun(@numel, acq_time)).';
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
boxplot(ax_hdl, bp_y, bp_x);
ax_hdl.XLabel.String = 'Layer';
ax_hdl.YLabel.String = 'Acquisition duration (s)';
%% Find tiles with very long acq time
vis_layer = [11, 12, 13, 14, 15];
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 3];

for i = 1 : numel(vis_layer)
    check_layer = vis_layer(i);
    check_tiles = tiles_cell{check_layer};
    acq_init = cat(1, check_tiles.t_init);
    [acq_init, acq_init_idx] = sortrows(acq_init);
    check_tiles_sorted = check_tiles(acq_init_idx);

    acq_duration_sorted = arrayfun(acq_time_fun, check_tiles_sorted);
    tmp_ax_hdl = subplot(5, 1, i);
    plot(tmp_ax_hdl, acq_duration_sorted);
    hold(tmp_ax_hdl, 'on');
    tmp_ax_hdl.Title.String = sprintf("Layer %d", check_layer);
    tmp_ax_hdl.YLabel.String = 'Acquisition duration (s)';
end
tmp_ax_hdl.XLabel.String = 'Ordered Acq Index';
fig_fp = fullfile(vis_folder, sprintf('Ordered_tile_acq_duration_11_to_15.png'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
% ax_hdl.YLabel.String = 'Acquisition duration (s)';
%%
test_tile = check_tiles_sorted(284);
test_smip = test_tile.load_step_mip();

seconds(datetime(test_tile.t_done) - datetime(test_tile.t_init))
%%
layer_idx = 15;
tv = WBIMTileViewer(check_tiles_sorted);
tv.update_tile_grid();
tv.init_step_mip(layer_data(layer_idx).smip{3}, layer_data(layer_idx).local_tile_info.step_mip_pixel_yxz_um);
tv.ui_init_ctrl();
% %%