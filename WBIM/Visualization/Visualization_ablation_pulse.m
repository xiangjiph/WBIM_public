DataManager = WBIMFileManager;
vis_folder = fullfile(DataManager.fp_processing_disk, 'Visualization', ...
    'Pipeline');
%% Tissue ablation 
stack_size = 200;
overlap_size = 100;
ablation_size = stack_size - overlap_size;
abl_skl_fine_lz_um = 40;
abl_skl_rough_dz_um = 5;
abl_skl_fine_dz_um = 3;
abl_skl_rough_abl_z_um = ablation_size - abl_skl_fine_lz_um;
abl_skl_fine_f_J_cm2 = 5;
abl_skl_rough_f_J_cm2 = 8;
%% Without amplification 
rough_abl_z = 0 : abl_skl_rough_dz_um : abl_skl_rough_abl_z_um;
rough_abl_f = repelem(abl_skl_rough_f_J_cm2, numel(rough_abl_z));
fine_abl_z = (abl_skl_rough_abl_z_um + abl_skl_rough_dz_um) : ...
    abl_skl_fine_dz_um : ablation_size;
fine_abl_f = repelem(abl_skl_fine_f_J_cm2, numel(fine_abl_z));
amp_1_l = 200;
amp_2_l = 100;
plt_x = [rough_abl_z, fine_abl_z];
plt_y = [rough_abl_f, fine_abl_f];
amp_coeff_1 = exp((ablation_size - plt_x) / amp_1_l);
amp_coeff_2 = exp((ablation_size - plt_x) / amp_2_l);
plt_y_1 = plt_y .* amp_coeff_1;
plt_y_1(1) = plt_y_1(1) * 1.5;
plt_y_2 = plt_y .* amp_coeff_2;
plt_y_2(1) = plt_y_2(1) * 1.5;

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
sc_hdl = scatter(ax_hdl, plt_x, plt_y, 'filled');
hold(ax_hdl, 'on');
sc_hdl_1 = scatter(ax_hdl, plt_x, plt_y_1, 'filled');
sc_hdl_2 = scatter(ax_hdl, plt_x, plt_y_2, 'filled');
% for i = 1 : numel(plt_x)
%     line(ax_hdl, plt_x([i, i]), [0, plt_y(i)], 'Color', 'r');
% end
ax_hdl.YLim(1) = 0;
view([0, 90]);
ax_hdl.XLabel.String = 'Z (\mum)';
ax_hdl.YLabel.String = 'Fluence (J/cm^2)';
legend(ax_hdl, [sc_hdl, sc_hdl_1, sc_hdl_2], 'After Scan', 'Refinement 1', 'Refinement 2', ...
    'Location', 'best');
ax_hdl.Title.String = 'Skull Ablation Scheme';

fig_fp = fullfile(vis_folder, 'Skull_ablation_scheme_3_rounds.png');
fun_print_image_in_several_formats(fig_hdl, fig_fp);

