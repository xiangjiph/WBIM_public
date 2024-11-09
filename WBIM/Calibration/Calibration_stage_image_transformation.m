exp_group = "Calibration";
exp_name = sprintf('WBIM%s', datestr(now, 'yyyymmdd'));
%% Tile Settings
imp = struct;
% xyz in the sample coordinate (fast, slow, z)
imp.scan.stack_xyz_size = [1024, 1280, 150];
imp.scan.pixel_xyz_size_um = [0.5, 0.5, 1];
imp.scan.overlap_xyz_size = [25, 25, 50];

imp.overview.stack_xyz_size = [700, 700, 1];
imp.overview.pixel_xyz_size_um = [1.0, 1.0, 0];
imp.overview.overlap_xyz_size = [0, 0, 0];
imp.active_channel = [1]; 
%% Ablation parameters
ablation_p = struct;
ablation_p.z_axis_step_um = 15; % To be determined by experiment 
ablation_p.site_scale_factor_fast = 1/sqrt(2);
ablation_p.site_scale_factor_slow = 1/2;
ablation_p.sample_mass_kg = 0.7;
% Need to offset the extra ablated thickness
ablation_p.ablation_thickness_um = (imp.scan.stack_xyz_size(3) - imp.scan.overlap_xyz_size(3)) * ...
    imp.scan.pixel_xyz_size_um(3);
ablation_p.max_acceleration_m_s2 = 0.5;
%% Initialization
if exist('wbim', 'var')
    wbim.delete();
    clearvars('wbim');
end
if exist('hSI', 'var')
    wbim = WBIMControl(exp_group, exp_name, imp, ablation_p, hSI);
else
    wbim = WBIMControl(exp_group, exp_name, imp, ablation_p, []);
end
%%
wbim.imaging.update_grid_parameters(imp);
wbim.init_tile_manager();
%%
wbim.add_current_pos_to_candidate_map();
wbim.tile_manager.expand_candidate_map(1, 'all');
wbim.set_microscope_mode(WBIMMicroscopeMode.Explore);
wbim.imaging.auto_exploration = false;
wbim.image_candidate_tiles();
%%
wbim.tile_manager.visualize_imaged_tile_mip_sample_yx_um(1);