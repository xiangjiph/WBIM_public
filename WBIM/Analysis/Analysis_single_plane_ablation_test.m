clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'TestSample';
exp_name = 'WBIM20230425001_20230511';
% im_root_folder = '/net/birdstore/Vessel/WBIM/Acquisition/TestSample/WBIM20230425001_20230511/00001/Ablation/Single_plane';
% abl_plane   = [nan, 0, 10, 15, 20, 25, 30, 35, 40, 45, 45, 50, 55, 60, 65, 65, 70, 75, 75, 80, 84];
% abl_fluence = [nan, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 10, 10, 15, 10, 10, 15, 10, 10];

im_root_folder = '/net/birdstore/Vessel/WBIM/Acquisition/TestSample/WBIM20230425001_20230511/00001/Ablation/Single_plane_site_2';
abl_plane   = [nan, 0:3:54];
abl_fluence = [nan, repelem(6, 1, 8), repelem(4.5, 1, 9), repelem(3, 1, 2)];

vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'single_plane_ablation_site_2');
%%
% All tiles are 700 x 700 x 200 um at 1 um isotropic resolution. 
tiff_list = dir(fullfile(im_root_folder, '*.tif'));
num_tiff = numel(tiff_list);
num_ch = 2;
im_stack_cell = cell(num_ch, num_tiff);
for i = 1 : num_tiff
    tmp_fp = fullfile(tiff_list(i).folder, tiff_list(i).name);
    tmp_im = DataManager.load_single_tiff(tmp_fp);
    for j = 1 : num_ch
        im_stack_cell{j, i} = medfilt3(tmp_im(:, :, j:num_ch:end));
    end    
end
im_stack_cell_zxy = cellfun(@(x) permute(uint16(x) * 2, [3, 2, 1]), im_stack_cell, 'UniformOutput', false);
%% Visualize the zy view 
vis_sec = 250;
vis_cell = cell(1, num_tiff);
for i = 1 : num_tiff
    tmp_stacks = im_stack_cell_zxy(:, i);
    tmp_im = cellfun(@(x) fun_stretch_contrast(x(:, :, vis_sec), 0.01 ,0.99), tmp_stacks, 'UniformOutput', false);
    tmp_vis_im = repmat(tmp_im{2}, 1, 1, 3) * 0.5;
    tmp_vis_im(:, :, 1) = max(tmp_vis_im(:, :, 1), tmp_im{1});        
    vis_cell{i} = tmp_vis_im;
end
%%
if ~isfolder(vis_folder)
    mkdir(vis_folder);
end
write_gif_Q = false;
fig_fp = fullfile(vis_folder, sprintf('%s_%s_zx_sec_%d.gif', exp_group, ...
    exp_name, vis_sec));
write_im_Q = true;
fig_hdl = figure;
t_hdl = tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');
ax_hdl = nexttile();
for i = 1 : num_tiff 
    if i == 1 
        im_hdl = imagesc(ax_hdl, vis_cell{i});
        hold(ax_hdl, 'on');
        ax_hdl.DataAspectRatio = [1,1,1];
        ax_hdl.XLabel.String = 'X (\mum)';
        ax_hdl.YLabel.String = 'Z (\mum)';
        ax_hdl.Title.String = 'Before ablation';
    else
        if i == 2
            l_hdl = plot(ax_hdl, [50, 650], abl_plane([i, i]), 'w-.', 'LineWidth', 1);
        else
            l_hdl.YData = abl_plane([i, i]);
        end
        im_hdl.CData = vis_cell{i};
        ax_hdl.Title.String = sprintf('%.1f J/cm^2 @ %d \\mum', ...
            abl_fluence(i), abl_plane(i));
    end    
    if write_gif_Q
        frame = getframe(fig_hdl);
        [imind, cm] = rgb2ind(frame2im(frame), 256);
        if i == 1
            imwrite(imind, cm, fig_fp, 'gif', 'LoopCount', inf);
        else
            imwrite(imind, cm, fig_fp, 'gif', 'WriteMode', 'append');
        end
    elseif write_im_Q
        fig_fp = fullfile(vis_folder, sprintf('zx_sec_%d', vis_sec), ...
            sprintf('%s_%s_zx_sec_%d_%02d.png', exp_group, exp_name, vis_sec, i));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);
    else
        pause(1);
    end
end
%%