% Generate checker board pattern 
grid_size = [7, 7];
grid_spacing_um = 250;
[sub1, sub2] = ndgrid(1 : grid_size(1), 1 : grid_size(2));
grid_ind = sub2ind(grid_size, sub1(:), sub2(:));
selected_Q = mod(grid_ind, 2) == 1;
grid_mask = false(grid_size);
grid_mask(grid_ind(selected_Q)) = true;
grid_mask_um = imresize(grid_mask, grid_spacing_um, 'nearest');
%%
roi_abl_ssf_fast = 4;
roi_abl_ssf_slow = 0.75;
roi_abl_z_step_um = 10;
bbox_z_l_um = 50;

roi_abl_peak_fluence_J_cm2 = 4;
current_xyz_um = wbim.sample_xyz_um;
current_fp_z_um = wbim.stage_fp_xyz_um(3) + ...
    (0 : roi_abl_z_step_um : bbox_z_l_um);

abl_path_2D = wbim.ablation.construct_path_for_sample_space_mask(...
    grid_mask_um, current_xyz_um([2,1]), current_fp_z_um, ...
    roi_abl_ssf_fast, roi_abl_ssf_slow, roi_abl_peak_fluence_J_cm2, false);
abl_path_2D.visualize_mask_with_trajectory();
%%
wbim.ablation.execute_ablation_paths(abl_path_2D);
%%
wbim.switch_to_imaging_mode();