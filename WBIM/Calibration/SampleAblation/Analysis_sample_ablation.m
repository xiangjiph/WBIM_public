DataManager = WBIMFileManager;
root_folder = fullfile(DataManager.DATA_ROOT_PATH, ...
    'Ablation_test\20220511');
analysis_folder = fullfile(root_folder, 'result');
pixel_size_um = [1, 1, 5];
im_fn = 'WBIM_abl_test_00005.tif';
[~, data_name, ~] = fileparts(im_fn);
vis_folder = fullfile(analysis_folder, data_name);
if ~isfolder(vis_folder)
    mkdir(vis_folder);
end
im_fp = fullfile(root_folder, im_fn);
im_hdl = ScanImageTiffReader.ScanImageTiffReader(im_fp);
im_data = im_hdl.data();
im_data = permute(im_data, [2,1,3]);
im_data = imresize3(im_data, size(im_data) .* pixel_size_um);
im_size = size(im_data);
im_info = im_hdl.metadata();
im_hdl.close();
%% Preprocessing
im_data = medfilt3(im_data);
im_data_zyx = permute(im_data, [3,1,2]);
%% Get the background connected component
max_bg_depth = 50;
est_bg_th = 5000;
bg_mask = im_data < est_bg_th;
bg_mask_sm = imclose(bg_mask, strel('sphere', 3));
bg_mask_cc = bwconncomp(bg_mask_sm);
[bg_mask_cc_size, cc_ind] = sort(cellfun(@numel, bg_mask_cc.PixelIdxList), 'descend');
bg_mask_sm(cat(1, bg_mask_cc.PixelIdxList{cc_ind(bg_mask_cc_size < 2e4)})) = false;
bg_mask_sm = imclose(bg_mask_sm, strel('cuboid', [25, 25, 2]));
% Find the boundary voxels in z direction
bg_mask_sm_g = (bg_mask_sm(:, :, 2:end) - bg_mask_sm(:, :, 1:end-1)) < 0;
% bg_mask_sm_g = padarray(bg_mask_sm_g, [0, 0, 1], 'post');
bg_mask_edge_ind = find(bg_mask_sm_g(:, :, 1:max_bg_depth));
bg_mask_edge_sub = fun_ind2sub(size(bg_mask_sm_g), bg_mask_edge_ind);

bg_mask_edge_sub_m = mean(bg_mask_edge_sub, 1);
% Assume z is a function of x and y
X_mat = cat(2, bg_mask_edge_sub(:, [1,2]), ones(numel(bg_mask_edge_ind), 1));
n_vec = inv(X_mat.' * X_mat) * (X_mat.' * bg_mask_edge_sub(:, 3));
bg_mask_sm_p = bg_mask_sm;
bg_mask_sm_p_ind = find(bg_mask_sm_p);
bg_mask_sm_p_sub = fun_ind2sub(im_size, bg_mask_sm_p_ind);
is_outlier_Q = (bg_mask_sm_p_sub(:, 3) - bg_mask_sm_p_sub(:, 1) .* n_vec(1) - ...
    bg_mask_sm_p_sub(:, 2) .* n_vec(2) - n_vec(3) ) > 2;
bg_mask_sm_p(bg_mask_sm_p_ind(is_outlier_Q)) = false;
%% Visualization
[plt_x, plt_y] = meshgrid(1:im_size(2), 1:im_size(1));
plt_z = plt_x .* n_vec(2) + plt_y .* n_vec(1) + n_vec(3);
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
surf(ax_hdl, plt_x, plt_y, plt_z);
hold(ax_hdl, 'on');
scatter3(ax_hdl, bg_mask_edge_sub(1:100:end, 2), bg_mask_edge_sub(1:100:end, 1), ...
    bg_mask_edge_sub(1:100:end, 3));
%%
x_list = 100 : 150 : 700;
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 2];
for i = 1 : numel(x_list)
    ax_hdl = subplot(3,2,i);
    imagesc(ax_hdl, im_data_zyx(:, :, x_list(i)));
    ax_hdl.DataAspectRatio = [1,1,1];
    ax_hdl.XLabel.String = 'Y (\mum)';
    ax_hdl.YLabel.String = 'Z (\mum)';
    ax_hdl.Title.String = sprintf('%d pulse(s) per site', i);
end
fig_fp = fullfile(vis_folder, sprintf('%s_ablation_site_im_zy.png', data_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

%% Determine the xy position of the ablation site:
% [662, 100, 1] * n_vec
abl_z = 50;
disp_dz = 20;
abl_p_im = im_data(:, :, abl_z);
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 2];
ax_hdl = subplot(2,2,1);
imagesc(ax_hdl, im_data(:, :, abl_z));
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.Title.String = 'Ablation plan image';
cbar_1 = colorbar(ax_hdl);
cbar_1.Label.String = 'Intensity';
ax_hdl_2 = subplot(2,2,2);
imagesc(ax_hdl_2, abl_z - max(0, plt_z));
cbar_2 = colorbar(ax_hdl_2);
cbar_2.Label.String = 'Depth(\mum)';
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_2.XLabel.String = 'X (\mum)';
ax_hdl_2.YLabel.String = 'Y (\mum)';
ax_hdl_2.Title.String = 'Distance to the bulk surface';

ax_hdl_3 = subplot(2,2,3);
imagesc(ax_hdl_3, im_data(:, :, abl_z - disp_dz));
ax_hdl_3.DataAspectRatio = [1,1,1];
ax_hdl_3.XLabel.String = 'X (\mum)';
ax_hdl_3.YLabel.String = 'Y (\mum)';
ax_hdl_3.Title.String = sprintf('Ablation plan image - %d \\mum', disp_dz);
cbar_1 = colorbar(ax_hdl_3);
cbar_1.Label.String = 'Intensity';

ax_hdl_4 = subplot(2,2,4);
imagesc(ax_hdl_4, im_data(:, :, abl_z + disp_dz));
ax_hdl_4.DataAspectRatio = [1,1,1];
ax_hdl_4.XLabel.String = 'X (\mum)';
ax_hdl_4.YLabel.String = 'Y (\mum)';
ax_hdl_4.Title.String = sprintf('Ablation plan image + %d \\mum', disp_dz);
cbar_1 = colorbar(ax_hdl_4);
cbar_1.Label.String = 'Intensity';
%%
fig_fp = fullfile(vis_folder, sprintf('%s_ablation_site_im_pm_20um.png', data_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);