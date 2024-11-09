%% Background fluctuation 
vis_sec = 645;
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
ax_1 = subplot(1,2,1);
imagesc(ax_1, smip_im{1}(:, :, vis_sec), 'CDataMapping', 'direct');
ax_1.ColorScale = 'linear';
ax_1.DataAspectRatio = [1,1,1];
c_bar_1 = colorbar(ax_1);
colormap(ax_1, 'jet');
ax_1.CLim(2) = 3e3;
ax_1.Title.String = sprintf('SHG section %d', vis_sec);
ax_1.YLabel.String = 'Y (5 \mum)';
ax_1.XLabel.String = 'X (5 \mum)';

ax_2 = subplot(1,2,2);
imagesc(ax_2, rescale(smip_im{1}(:, :, vis_sec)).^(1/3));
ax_2.DataAspectRatio = [1,1,1];
ax_2.ColorScale = 'linear';
c_bar = colorbar(ax_2);
% c_bar.Limits(1) = 10;
% c_bar.Limits(2) = 1e4;
colormap(ax_2, 'jet');
ax_2.Title.String = sprintf('SHG section %d, enhanced \\gamma = 1/3', vis_sec);
c_bar.Label.String = 'Normalized Int';
fig_fp = fullfile(vis_folder, sprintf('SMIP_CH1_sec_%d_gamma_1_3_background_fluctuation.png',...
    vis_sec));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Check missing tiles
% missing data started at 638
missing_sec = 638;
channel_list = tile_str{1}(1).channel;
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 3;
for i = 1 : 3
    ax_hdl = subplot(1,3,i);
    imagesc(smip_im{i}(:, :, missing_sec));
    ax_hdl.DataAspectRatio = [1,1,1];
    ax_hdl.CLim(2) = 1.5e4;
    ax_hdl.Title.String = sprintf('CH %d', channel_list(i));
end
fig_fp = fullfile(vis_folder, sprintf('SMIP_sec_%d_missing_ch4.png', vis_sec));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% The last tile with signal in CH4
layer_idx = 14;
tile_idx = 770;
ch_idx = 3;
last_tile = tile_str{layer_idx}(tile_idx);
last_mip = last_tile.load_step_mip();
last_mip = last_mip{ch_idx};
last_mip_zyx = permute(last_mip, [3, 1, 2]);

fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
ax_hdl = subplot(1,2,1);
imagesc(max(last_mip_zyx, [], 3));
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.Title.String = sprintf('Tile %d', tile_idx);
ax_hdl.CLim(2) = 4e3;
colormap(ax_hdl, 'jet');

last_stat = last_tile.load_stat();
last_stat = last_stat.CH3;
ax_hdl_2 = subplot(1,2,2);
plot(ax_hdl_2, last_stat.z_stat.mean);
ax_hdl_2.YLabel.String = 'Average Z intensity';
ax_hdl_2.XLabel.String = 'Z (\mum)';
fig_fp = fullfile(vis_folder, sprintf('SMIP_CH3_layer_%d_tile_%d_z_int.png', layer_idx, tile_idx));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Merge two channels
test_sec = 100;
im_1 = smip_im{2}(:, :, test_sec);
im_2 = smip_im{3}(:, :, test_sec);
% There are still saturated regions in im_2... 
% Saturated in im 1
im_1_min_int = 0;
im_sat_int = 30000;
m_s_1 = im_1 > im_sat_int;
m_m_1 = im_1 < im_sat_int & im_1 > im_1_min_int;
m_ns_2 = im_2 < im_sat_int;
% Linear fit between the intensity of the median-intensity voxels in
% channel 1
int_m_1 = double(im_1(m_m_1));
int_m_2 = double(im_2(m_m_1));

bin_int_edge = 0 : 100 : im_sat_int;
[bin_ind, bin_val] = fun_bin_data_to_idx_list_by_edges(int_m_2, bin_int_edge);
bin_stat = fun_analysis_get_basic_statistics_in_bins(int_m_1, bin_ind);
is_valid_Q = arrayfun(@(x) ~isempty(x.mean), bin_stat); 
int_m_1_q3 = arrayfun(@(x) x.prctile_val([7, 8, 9]), bin_stat(is_valid_Q), 'UniformOutput', false);
int_m_1_q3 = cat(1, int_m_1_q3{:});
int_m_1_e1 = int_m_1_q3(:, 2) - int_m_1_q3(:, 1);
int_m_1_e2 = int_m_1_q3(:, 3) - int_m_1_q3(:, 2);
int_m_1_aen = (int_m_1_e1 + int_m_1_e2) ./ int_m_1_q3(:, 2);

bin_val_ne = bin_val(is_valid_Q);
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_1 = axes(fig_hdl);
histogram2(int_m_2, int_m_1, 'DisplayStyle', 'tile');
ax_1.DataAspectRatio = [1,1,1];
ax_1.XLim = [0, im_sat_int];
ax_1.YLim = [im_1_min_int, im_sat_int];
ax_1.ColorScale = 'log';
ax_1.XLabel.String = 'CH4';
ax_1.YLabel.String = 'CH3';
ax_1.Title.String = sprintf('%s_%s_SMIP_sec_%d', exp_group, exp_name, test_sec);
ax_1.Title.Interpreter = 'none';
c_bar = colorbar(ax_1);
c_bar.Label.String = 'Counts';
hold(ax_1, 'on');

errorbar(ax_1, bin_val_ne, int_m_1_q3(:, 2), int_m_1_q3(:, 2) - int_m_1_q3(:, 1), ...
    int_m_1_q3(:, 3) - int_m_1_q3(:, 2), 'LineWidth', 0.5);
% Somehow the width is broden 
%%
cut_off_int = 10e3;
selected_bin_Q = bin_val_ne < cut_off_int;
selected_pxl_Q = int_m_2 < cut_off_int;
fit_direct = fitlm(int_m_2(selected_pxl_Q), int_m_1(selected_pxl_Q));
fit_med = fitlm(bin_val_ne(selected_bin_Q), int_m_1_q3(selected_bin_Q, 2));
fprintf('Cut off intensity: %d\n', cut_off_int);
disp(fit_med);
plt_x = linspace(ax_1.XLim(1), cut_off_int, 1e3);
plt_y = plt_x * fit_med.Coefficients.Estimate(2) + fit_med.Coefficients.Estimate(1);
fit_line_hdl = plot(ax_1, plt_x, plt_y, 'LineWidth', 1.5);
leg_txt = sprintf('Fit to median up to %.2e\nk = %.3e \\pm %.2e\nb = %.3e \\pm %.2e\nR^2-A = %.3f', ...
    cut_off_int, ...
    fit_med.Coefficients.Estimate(2), fit_med.Coefficients.SE(2), ...
    fit_med.Coefficients.Estimate(1), fit_med.Coefficients.SE(1), ...
    fit_med.Rsquared.Adjusted);
legend(ax_1, fit_line_hdl, leg_txt, 'Location', 'southeast');

% fig_fp = fullfile(vis_folder, sprintf('%s_%s_smip_sec_%d_ch_3_vs_4_fit_up_to_%d.png', ...
%     exp_group, exp_name, test_sec, cut_off_int));
% fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
int_scale_k = 2.5416;
int_scale_b = 23;
% saturated_mask = im_1 > 2.5e4;
im_2_e = double(im_2) * int_scale_k + int_scale_b;
saturated_mask = im_2_e > im_sat_int;
im_c = double(im_1);
im_c(saturated_mask) = max(im_c(saturated_mask), im_2_e(saturated_mask));
im_c(~saturated_mask) = (im_c(~saturated_mask) + im_2_e(~saturated_mask)) ./ 2;
im_c = im2uint16(rescale(im_c));

fig_hdl = figure;
ax_1 = subplot(1,3,1);
imagesc(ax_1, im_1);
ax_2 = subplot(1,3,2);
imagesc(ax_2, im_2);
ax_3 = subplot(1,3,3);
imagesc(ax_3, im_c);
[ax_1.DataAspectRatio, ax_2.DataAspectRatio, ax_3.DataAspectRatio] = deal([1,1,1]);
linkaxes([ax_1, ax_2, ax_3], 'xy');

% Merge two channels

%%
% figure;imshowpair(im_1, im_c);
figure;
histogram(im_c);
set(gca, 'YScale', 'log');







































