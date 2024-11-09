exp_group = "Dryrun";
exp_name = sprintf('Test_%s', datestr(now, 'yyyymmdd'));
%% Tile Settings
imp = struct;
% xyz in the sample coordinate (fast, slow, z)
% Align the last plane of the scan mode and explore mode 
img_wavelength_nm = 950;
img_gdd_fs2 = 15000;
scan_para = struct;
scan_para.z_offset_um = 0;
scan_para.stack_xyz_size = [896, 1536, 100];
scan_para.pixel_xyz_size_um = [0.25, 0.25, 1];
scan_para.stack_xyz_size_um = scan_para.stack_xyz_size .* scan_para.pixel_xyz_size_um;
scan_para.overlap_xyz_size = [64, 64, 60];
scan_para.ablation_thickness_um = (scan_para.stack_xyz_size(3) - ...
    scan_para.overlap_xyz_size(3)) * scan_para.pixel_xyz_size_um(3);
scan_para.imaging_power = {'none', 0.08};

exp_para = struct;
exp_para.pixel_xyz_size_um = [1.0, 1.0, 5];
exp_para.z_offset_um = - scan_para.ablation_thickness_um;
exp_para.image_thickness_um = (scan_para.stack_xyz_size_um(3) + scan_para.ablation_thickness_um + ...
    exp_para.pixel_xyz_size_um(3));
exp_para.stack_xyz_size_um = [700, 700, exp_para.image_thickness_um];
exp_para.stack_xyz_size = round(exp_para.stack_xyz_size_um ./ exp_para.pixel_xyz_size_um);
exp_para.overlap_xyz_size = [0, 0, 0];
exp_para.imaging_power  = {'none', 0.10};

imp.scan = scan_para;
imp.explore = exp_para;
imp.active_channel = [1, 2]; 
imp.bg_detection_channel = [1, 2];
imp.im_qc_channel = [1, 2];
%% Ablation parameters
target_ablation_thickness = (imp.scan.stack_xyz_size(3) - imp.scan.overlap_xyz_size(3)) * ...
    imp.scan.pixel_xyz_size_um(3);

abl_detection_ch = [1, 2];
with_skull_Q = true;
% Without diffuser
abl_p = {};
abl_p{end+1} = WBIMSPAblation(...
    'material', [WBIMAblationMaterial.Tissue, WBIMAblationMaterial.Vessel, ...
    WBIMAblationMaterial.Membrane, WBIMAblationMaterial.TissueInSkull], ...
    'with_diffuser_Q', false, 'site_scale_factor_fast', 2, 'site_scale_factor_slow', 0.75, ...
    'z_axis_step_um', 5, 'peak_fluence_J_cm2', 1.5, 'abl_z_offset_r_um', 0, ...
    'abl_range_r_um', [0, target_ablation_thickness]);

abl_p{end+1} = WBIMSPAblation(...
    'material', [WBIMAblationMaterial.Bone], ...
    'with_diffuser_Q', false, 'site_scale_factor_fast', 1, 'site_scale_factor_slow', 0.5, ...
    'z_axis_step_um', 2, 'peak_fluence_J_cm2', 1, 'abl_z_offset_r_um', 0, ...
    'abl_range_r_um', [0, target_ablation_thickness]);

abl_p = cat(1, abl_p{:});
%% Initialization
if exist('wbim', 'var')
    wbim.delete();
    clearvars('wbim');
end
if exist('hSI', 'var')
    wbim = WBIMControl(exp_group, exp_name, imp, abl_p, hSI);
else
    wbim = WBIMControl(exp_group, exp_name, imp, abl_p, []);
end
%%
% wbim.imaging.update_grid_parameters(imp);