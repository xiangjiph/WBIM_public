%% Compare sample preparation methods
DataManager = WBIMFileManager;
exp_group_name = 'TestSample';
exp_name = 'Clearing_test';
root_data_folder = DataManager.fp_experiment(exp_group_name, exp_name);
%% DISCO-PBS
% The second channel is DAPI
% The 525/40 channel contains autofluorescence signal 
set_folder_name = 'ZW001_DAPI_PBS';
set_folder_path = fullfile(root_data_folder, set_folder_name);
im_fn = 'WBIM_00000.tif';

im_fp = fullfile(set_folder_path, im_fn);
im_disco_pbs = uint16(DataManager.load_single_tiff(im_fp)) * 2;
im_disco_pbs_2 = im_disco_pbs(:, :, 2:2:end);

im_disco_pbs_2_zyx = permute(im_disco_pbs_2, [3, 1, 2]);
%% DISCO-ECi
set_folder_path = fullfile(root_data_folder, 'ZW001_DAPI_ECi');
im_fn = 'ECi_00003.tif';
im_fp = fullfile(set_folder_path, im_fn);
im_disco_eci = uint16(DataManager.load_single_tiff(im_fp)) * 2;
im_disco_eci_zyx = permute(im_disco_eci, [3,1,2]);
%% Brain without the skull
set_folder_name = 'DK98_47TDE';
set_folder_path = fullfile(root_data_folder, set_folder_name);
im_fn = 'WBIM_00010.tif';
im_fp = fullfile(set_folder_path, im_fn);
im_tde_vsl = uint16(DataManager.load_single_tiff(im_fp)) * 2;

%% 47% TDE with skull
% Image through the skull
set_folder_name = '20220627001_47TDE';
set_folder_path = fullfile(root_data_folder, set_folder_name);
im_fn = 'WBIM_00006.tif';
im_fp = fullfile(set_folder_path, im_fn);
im_tde_skl = uint16(DataManager.load_single_tiff(im_fp)) * 2;
im_tde_skl_vsl = im_tde_skl(:, :, 1:2:end);
im_tde_skl_shg = im_tde_skl(:, :, 2:2:end);
%% 60% DMSO 40% OptiPrep
set_folder_name = 'WBIM20220629002_20230320';
set_folder_path = fullfile(root_data_folder, set_folder_name);
im_fn = 'WBIM_00000.tif';
im_fp = fullfile(set_folder_path, im_fn);
im_d6o4 = uint16(DataManager.load_single_tiff(im_fp)) * 2;
im_d6o4_vsl = im_d6o4(:, :, 1:2:end);
im_d6o4_shg = im_d6o4(:, :, 2:2:end);
%% 100% TDE 1.0 na
set_folder_name = 'WBIM20220623001_20230320';
set_folder_path = fullfile(root_data_folder, set_folder_name);
im_fn = 'WBIM_00004.tif';
im_fp = fullfile(set_folder_path, im_fn);
im_tde100 = uint16(DataManager.load_single_tiff(im_fp)) * 2;
im_tde100_vsl = im_tde100(:, :, 1:2:end);
im_tde100_shg = im_tde100(:, :, 2:2:end);
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
imagesc(ax_hdl, max(medfilt3(im_tde100_vsl(:, :, 1 : 300)), [], 3));
ax_hdl.DataAspectRatio = [1,1,1];
colormap("gray");
ax_hdl.XAxis.Visible = 'off';
ax_hdl.YAxis.Visible = 'off';
vis_title = "100% TDE Vessel";
[~, im_fn_woe, ~] = fileparts(im_fn);
fp_pf = strrep(vis_title, ' ', '_');
fig_fp = fullfile(set_folder_path, 'vis', sprintf('%s_%s_%s_300um_mip.png', set_folder_name, im_fn_woe, fp_pf));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
% vis_stack = im_disco_pbs_2;
% vis_z = [25, 75, 150, 250];

% vis_stack = im_disco_eci;
% vis_z = [25, 200, 400, 800];

% vis_stack = im_tde_vsl;
% vis_z = [100, 250, 400, 600];
% vis_z = [100, 200, 300, 400];
% vis_title = "047TDE Vessel";

% vis_stack = im_tde_skl_shg;
% vis_z = [50, 100, 150, 200];

% vis_stack = im_tde_skl_shg;
% vis_z = [50, 150, 200, 250];
% vis_title = "047TDE Skull";

% vis_stack = im_tde_skl_vsl;
% vis_z = [50, 150, 200, 250];
% vis_title = "047TDE Vessel";

% vis_stack = im_d6o4_vsl;
% vis_z = [100, 200, 250, 300];
% vis_title = "60% DMSO 40% Optiprep Vessel";

% vis_stack = im_d6o4_shg;
% vis_z = [50, 100, 150, 200];
% vis_title = "60% DMSO 40% Optiprep SHG";

vis_stack = im_tde100_vsl;
vis_z = [100, 200, 300, 400];
vis_title = "100% TDE Vessel";

[~, im_fn_woe, ~] = fileparts(im_fn);
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) .* 4;
num_im = numel(vis_z);
for i = 1 : num_im
    tmp_ax = nexttile();
    imagesc(tmp_ax, vis_stack(:, :, vis_z(i)));
    tmp_ax.DataAspectRatio = [1,1,1];
    tmp_ax.Title.String = sprintf('z = %d \\mum', vis_z(i));
    colormap('gray')
    tmp_ax.XAxis.Visible = 'off';
    tmp_ax.YAxis.Visible = 'off';
    tmp_ax.CLim = [0, 2^16-1];
end

fp_pf = strrep(vis_title, ' ', '_');
fig_fp = fullfile(set_folder_path, 'vis', sprintf('%s_%s_%s.png', set_folder_name, im_fn_woe, fp_pf));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% ZY Plane 
% vis_im = im_disco_pbs_2_zyx(:, :, 350);
% im_title = 'IDISCO PBS DAPI';
vis_im = im_disco_eci_zyx(:, :, 350);
im_title = 'IDISCO ECi DAPI';


[~, im_fn_woe, ~] = fileparts(im_fn);
fig_fp = fullfile(set_folder_path, 'vis', sprintf('%s_%s_zy.png', set_folder_name, im_fn_woe));

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
imagesc(ax_hdl, vis_im);
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Z (\mum)';
cb_hdl = colorbar(ax_hdl);
ax_hdl.CLim(1) = 0;
ax_hdl.Title.String = im_title;

fun_print_image_in_several_formats(fig_hdl, fig_fp);