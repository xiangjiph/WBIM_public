clc;clear;close all;
%%
DataManager = WBIMFileManager;
exp_group = 'TestSample';
exp_name = 'WBIM20220616002_20221027';
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'abl_roi_detection');
%%
tile_str = DataManager.load_tile_in_experiment(exp_group, exp_name);
%%
layer_idx = 8;
exp_tile = tile_str.Explore{layer_idx};
scan_tile = tile_str.Scan{layer_idx};

process_tile = scan_tile;
[smip_cell, local_tile_info] = WBIMTileManager.get_stitched_step_mip(...
    process_tile, [1,2]);
%% Visualize selected layer - raw data
vis_sec = 16;
fig_prefix = sprintf('%s_%s_layer_%d_%s_sec_%d', exp_group, exp_name, layer_idx, ...
    'Scan', vis_sec);
vis_shg = smip_cell{2}(:, :, vis_sec);
vis_vsl = smip_cell{1}(:, :, vis_sec);
vis_merge = fun_merge_image_stacks({uint16(vis_vsl), uint16(vis_shg)}, 'method', 'VesselBone');
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [4, 2];
a1 = nexttile;
imagesc(a1, vis_shg);
a1.DataAspectRatio = [1,1,1];
a1.XLabel.String = 'X/(5\mum)';
a1.YLabel.String = 'Y/(5\mum)';
colormap(a1, 'gray');
a1.Title.String = 'SHG';
a2 = nexttile;
imagesc(a2, vis_vsl);
a2.DataAspectRatio = [1,1,1];
a2.XLabel.String = 'X/(5\mum)';
a2.YLabel.String = 'Y/(5\mum)';
colormap(a2, 'gray');
a2.Title.String = 'Vessel';
a3 = nexttile;
image(a3, vis_merge);
a3.DataAspectRatio = [1,1,1];
a3.XLabel.String = 'X/(5\mum)';
a3.YLabel.String = 'Y/(5\mum)';
a3.Title.String = 'Merged';
fig_fp = fullfile(vis_folder, sprintf('%s_smip.png', fig_prefix));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Visualize selected layer - spectral unmixing
% run after obj.compute_single_channel_masks - smio_spectral_unmixing
layer_idx = 8;
vis_sec = 8;
fig_prefix = sprintf('%s_%s_layer_%d_%s_sec_%d', exp_group, exp_name, layer_idx, ...
    'Scan', vis_sec);

vis_shg = obj.smip_cell{2}(:, :, vis_sec);
vis_vsl = obj.smip_cell{1}(:, :, vis_sec);
vis_merge = fun_merge_image_stacks({uint16(vis_vsl), uint16(vis_shg)}, 'method', 'VesselBone');
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [4, 2];
a1 = nexttile;
imagesc(a1, vis_shg);
a1.DataAspectRatio = [1,1,1];
a1.XLabel.String = 'X/(5\mum)';
a1.YLabel.String = 'Y/(5\mum)';
colormap(a1, 'gray');
a1.Title.String = 'SHG';
a2 = nexttile;
imagesc(a2, vis_vsl);
a2.DataAspectRatio = [1,1,1];
a2.XLabel.String = 'X/(5\mum)';
a2.YLabel.String = 'Y/(5\mum)';
colormap(a2, 'gray');
a2.Title.String = 'Vessel';
a3 = nexttile;
image(a3, vis_merge);
a3.DataAspectRatio = [1,1,1];
a3.XLabel.String = 'X/(5\mum)';
a3.YLabel.String = 'Y/(5\mum)';
a3.Title.String = 'Merged';
fig_fp = fullfile(vis_folder, sprintf('%s_smip_sumix.png', fig_prefix));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% 
abl_vol_dect = WBIMAVD(smip_cell, local_tile_info, 'acquisition_mode', ...
    WBIMMicroscopeMode.Scan, 'has_skull_Q', true);
abl_vol_dect.compute_single_channel_masks();
abl_vol_dect.construct_labeled_mask('detect_vsl_in_skull_Q', false);
%% Visualize mask 
vis_vsl_m = abl_vol_dect.mask_cell{1}(:, :, vis_sec);
vis_shg_m = abl_vol_dect.mask_cell{2}(:, :, vis_sec);
vis_label_m = abl_vol_dect.labeled_mask(:, :, vis_sec);
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [4, 2];
a1 = subplot(1,3,1);
imshowpair(vis_vsl, vis_vsl_m);
a1.XAxis.Visible = 'on';
a1.YAxis.Visible = 'on';
a1.XLabel.String = 'X/5\mum';
a1.YLabel.String = 'Y/5\mum';
a1.Title.String = 'Vascular channel overlaid with mask';
a2 = subplot(1,3,2);
imshowpair(vis_shg, vis_shg_m);
a2.Title.String = 'SHG channel overlaid with mask';
a3 = subplot(1,3,3);
imagesc(a3, vis_label_m);
a3.DataAspectRatio = [1,1,1];
a3.Title.String = 'Labeled array';
a3.XLabel.String = 'X/5\mum';
a3.YLabel.String = 'Y/5\mum';
colormap(a3, 'jet');
fig_fp = fullfile(vis_folder, sprintf('%s_masks.png', fig_prefix));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Ablation parameters
%% Tile Settings
imp = struct;
% xyz in the sample coordinate (fast, slow, z)
% Align the last plane of the scan mode and explore mode 
% Also need to set the imaging wavelength 
scan_para = struct;
scan_para.z_offset_um = 0;
scan_para.stack_xyz_size = [896, 1536, 100];
scan_para.pixel_xyz_size_um = [0.25, 0.25, 1];
scan_para.stack_xyz_size_um = scan_para.stack_xyz_size .* scan_para.pixel_xyz_size_um;
scan_para.overlap_xyz_size = [64, 64, 0];
scan_para.ablation_thickness_um = (scan_para.stack_xyz_size(3) - ...
    scan_para.overlap_xyz_size(3)) * scan_para.pixel_xyz_size_um(3);
scan_para.imaging_power = {'none', 0.07};

imp.scan = scan_para;
% Ablation parameters
target_ablation_thickness = (imp.scan.stack_xyz_size(3) - imp.scan.overlap_xyz_size(3)) * ...
    imp.scan.pixel_xyz_size_um(3);

abl_vsl_fine_lz_um = 40;
abl_vsl_rough_dz_um = 5;
abl_vsl_fine_dz_um = 5;
abl_vsl_rough_abl_z_um = target_ablation_thickness - abl_vsl_fine_lz_um;

abl_skl_fine_lz_um = 40;
abl_skl_rough_dz_um = 5;
abl_skl_fine_dz_um = 5;
abl_skl_rough_abl_z_um = target_ablation_thickness - abl_skl_fine_lz_um;

abl_detection_ch = [1,2];
with_skull_Q = true;
% Without diffuser
abl_p = {};
abl_p{end+1} = WBIMSPAblation(...
    'material', [WBIMAblationMaterial.Tissue, WBIMAblationMaterial.Vessel, WBIMAblationMaterial.Membrane], ...
    'with_diffuser_Q', false, 'site_scale_factor_fast', 4, 'site_scale_factor_slow', 0.75, ...
    'z_axis_step_um', abl_vsl_rough_dz_um, 'peak_fluence_J_cm2', 5, 'abl_z_offset_r_um', 0, ...
    'abl_range_r_um', [0, abl_vsl_rough_abl_z_um]);

abl_p{end+1} = WBIMSPAblation(...
    'material', [WBIMAblationMaterial.Tissue, WBIMAblationMaterial.Vessel, WBIMAblationMaterial.Membrane], ...
    'with_diffuser_Q', false, 'site_scale_factor_fast', 4, 'site_scale_factor_slow', 0.75, ...
    'z_axis_step_um', abl_vsl_fine_dz_um, 'peak_fluence_J_cm2', 5, 'abl_z_offset_r_um', 0, ...
    'abl_range_r_um', [abl_vsl_rough_abl_z_um + abl_vsl_rough_dz_um, target_ablation_thickness]);

abl_p{end+1} = WBIMSPAblation(...
    'material', [WBIMAblationMaterial.TissueInSkull], ...
    'with_diffuser_Q', false, 'site_scale_factor_fast', 4, 'site_scale_factor_slow', 0.75, ...
    'z_axis_step_um', abl_vsl_fine_dz_um, 'peak_fluence_J_cm2', 3, 'abl_z_offset_r_um', 0, ...
    'abl_range_r_um', [0, target_ablation_thickness]);

abl_p{end+1} = WBIMSPAblation(...
    'material', [WBIMAblationMaterial.Bone], ...
    'with_diffuser_Q', false, 'site_scale_factor_fast', 1, 'site_scale_factor_slow', 0.5, ...
    'z_axis_step_um', abl_skl_rough_dz_um, 'peak_fluence_J_cm2', 10, 'abl_z_offset_r_um', 0, ...
    'abl_range_r_um', [0, abl_skl_rough_abl_z_um]);

abl_p{end+1} = WBIMSPAblation(...
    'material', [WBIMAblationMaterial.Bone], ...
    'with_diffuser_Q', false, 'site_scale_factor_fast', 1, 'site_scale_factor_slow', 0.5, ...
    'z_axis_step_um', abl_skl_fine_dz_um, 'peak_fluence_J_cm2', 10, 'abl_z_offset_r_um', 0, ...
    'abl_range_r_um', [abl_skl_rough_abl_z_um + abl_skl_rough_dz_um, target_ablation_thickness]);

abl_p = cat(1, abl_p{:});
%%
abl_vol_dect.construct_mask_yx_um_itp();
abl_vol_dect.construct_merged_ablation_mask(abl_p);
%%
wbim_config = WBIMConfig();
num_abl_paras = numel(abl_vol_dect.ablation_parameter);
abl_str = cell(num_abl_paras, 1);
for i = 1 : num_abl_paras
    tmp_abl_para = abl_vol_dect.ablation_parameter(i);
    mask_yx_um = abl_vol_dect.ablation_mask{i};
    tmp_abl_z_r_um = abl_vol_dect.ablation_z_r_um{i};
    tmp_abl_fp_z_um = abl_vol_dect.abs_fp_z_um_mm + tmp_abl_z_r_um;
    mask_center_yx_um = abl_vol_dect.region_bbox_ctr_yx_um;
    
    
    mask_xy_um = permute(wbim_config.sample_to_stage_im_yx_um(mask_yx_um), [2,1,3]);
    size_xy_um = size(mask_xy_um, [1,2]);
    % Pad 0 in z - in our setup we assume the sample to stage
    % transformation only changes xy coordinates
    center_stage_xy_um = wbim_config.sample_to_stage_xyz_um([mask_center_yx_um(2), ...
        mask_center_yx_um(1), 0]).';
    center_stage_xy_um = center_stage_xy_um(1:2);
    mask_xy_mm_um = round(center_stage_xy_um - size_xy_um/2 + 1);
    mask_xy_xx_um = mask_xy_mm_um + size_xy_um - 1;
    stage_bbox_xy_um = [mask_xy_mm_um, mask_xy_xx_um];

    abl_para = WBIMSPAblation.compute_ablation_parameters_for_multiple_planes(...
        tmp_abl_para.site_scale_factor_fast, tmp_abl_para.site_scale_factor_slow,...
        tmp_abl_para.with_diffuser_Q, tmp_abl_para.peak_fluence_J_cm2);
    
    abl_vol_3D = WBIMAblationPath3D(abl_para.fast_axis_speed_um_s, ...
        abl_para.fast_acceleration_length_um, abl_para.slow_axis_step_um, ...
        tmp_abl_fp_z_um, abl_para.power_W, abl_para.diffuserQ, stage_bbox_xy_um, mask_xy_um);
    abl_vol_3D.construct_path_from_mask();
    abl_path_2D = abl_vol_3D.path;
    abl_str{i} = abl_path_2D;
end
% abl_str = cat(1, abl_str{:});
for abl_i = 1 : num_abl_paras
    tmp_hdl = abl_str{abl_i}(vis_sec).visualize_mask_with_trajectory();
    tmp_fp = fullfile(vis_folder, sprintf('%s_abl_path_ch_%d.png', fig_prefix, abl_i));
    fun_print_image_in_several_formats(tmp_hdl, tmp_fp);
end
%%
abl_str_1 = cat(2, abl_str{:});
path_xyz = arrayfun(@(x) cat(2, x.trajectory(:, 1:2), repelem(x.abl_fp_abs_z_um, ...
    size(x.trajectory, 1), 1)), abl_str_1, 'UniformOutput', false);

path_xyz_1 = cat(1, path_xyz{15:17, 1});
path_xyz_2 = cat(1, path_xyz{15:17, 2});
min_z_um = min(path_xyz_1(:, 3));
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 2;
ax_hdl = axes(fig_hdl);
plot3(ax_hdl, path_xyz_1(:, 1), path_xyz_1(:, 2), path_xyz_1(:, 3) - min_z_um);
hold(ax_hdl, 'on');
plot3(ax_hdl, path_xyz_2(:, 1), path_xyz_2(:, 2), path_xyz_2(:, 3) - min_z_um);
ax_hdl.ZDir = 'reverse';
% ax_hdl.ZLim(1) = ax_hdl.ZLim(1) - 1;
legend(ax_hdl, 'Tissue', 'Skull');
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.DataAspectRatio = [400, 400, 1];

fig_fp = fullfile(vis_folder, sprintf('%s_abl_path_3D_15_to_17.png', fig_prefix));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
% abl_vol_dect.ablation_path = obj.construct_ablation_path_from_AVD(...
%     abl_vol_dect);
%% Load tiles and compute masks
test_data_fp = './Testing/Test_data/Ablation_path_planning.mat';
test_data = load(test_data_fp);
abl_para = test_data.abl_para;
abl_fp_z_abs_um = test_data.abl_fp_z_abs_um;
peak_fluence_J_cm2 = test_data.peak_fluence_J_cm2;
mask_xy_um = test_data.mask_xy_um;
stage_bbox_xy_um = test_data.stage_bbox_xy_um;
%%
i = 8;
abl_path_2D = WBIMAblationPath2Dv2(abl_para.fast_axis_speed_um_s, ...
    abl_para.fast_acceleration_length_um, abl_para.slow_axis_step_um, ...
    abl_fp_z_abs_um(i), peak_fluence_J_cm2(i), stage_bbox_xy_um, ...
    mask_xy_um(:, :, i));










