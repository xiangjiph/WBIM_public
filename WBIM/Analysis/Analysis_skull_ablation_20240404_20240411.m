clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'TestSample';
exp_name = 'WBIM20240305001_20240410';
data_root_folder = fullfile(DataManager.fp_layer(exp_group, exp_name, 20), ...
    'Explore');
si_im_cvrt = @(x) 2 * uint16(x);

% /nfs/birdstore/Vessel/WBIM/Acquisition/TestSample/WBIM20240305001_20240410/00020/Explore/Site_8

vis_root_folder = fullfile(DataManager.fp_visualization_folder(exp_group, exp_name), ...
    'Continuous_ablation');
vis_file_prefix = sprintf('%s_%s_cont_abl', exp_group, exp_name);
%% Iterative imaging and ablation
site_idx = 5;
site_folder = sprintf('Site_%d', site_idx);
exp_folder = fullfile(data_root_folder, site_folder);
exp_files = dir(fullfile(exp_folder, '*.tif'));

exp_fig_folder = fullfile(vis_root_folder, site_folder);
exp_fig_prefix = sprintf('%s_cont_abl_%s', vis_file_prefix, site_folder);

file_name_list = {exp_files.name};
file_count_list = cellfun(@(x) str2double(x(6:10)), file_name_list);
assert(issorted(file_count_list, 'ascend'));
num_files = numel(file_count_list);
% Load all files 
exp_im = arrayfun(@(x) si_im_cvrt(DataManager.load_single_tiff(fullfile(exp_folder, ...
    x.name))), exp_files, 'UniformOutput', false);

im_size = size(exp_im{1});
num_ch = 2;
pixel_size = 1;

registration_Q = true;
%%
% Smooth the image using median filter along the y axis
% exp_im_sm = cellfun(@(x) medfilt3(x, [3, 3, 1]), exp_im, 'UniformOutput', false);
fun_process_single_stack = @(y) cellfun(@(x) medfilt3(x, [3,3,3]), fun_im_deinterleave(y, 2), 'UniformOutput', false).';

exp_im_avg = cellfun(fun_process_single_stack, exp_im, 'UniformOutput', false);
exp_im_avg = cat(1, exp_im_avg{:});
exp_im_avg = cellfun(@(x) x(:, :, 1:end-1), exp_im_avg, 'UniformOutput', false);
% Move the tile after ablation up? 
%% Registration - for potential global deformation 
% Global translation is relatively small (~ 1 um in y and z). Neglect them
% at the moment.
reg_ch = 1;
non_reg_ch = 2;
exp_im_reg = exp_im_avg;
num_sets = size(exp_im_avg, 1);
disp_vec = zeros(3, num_sets, num_ch);
reg_vis_sec = 350;
if registration_Q
    mfft_th = 1e4;
    fft_search_range = [1, 1, 15];
    ref_im = exp_im_reg{1, reg_ch};
    for i = 2 : num_sets
        tmp_reg = exp_im_reg{i, reg_ch};
        tmp_ref_mask = ref_im >= mfft_th;
        tmp_reg_mask = tmp_reg >= mfft_th;
        [tmp_vec, tmp2] = MaskedTranslationRegistration(ref_im, tmp_reg, ...
            tmp_ref_mask, tmp_reg_mask, fft_search_range);
        % Remove translation registration? 
        tmp_vec(1:2) = 0;
        if any(tmp_vec)
            tmp_reg_mv = imtranslate(tmp_reg, tmp_vec);
            tmp_empty_mask = (tmp_reg_mv == 0);
            tmp_reg_mv(tmp_empty_mask) = tmp_reg(tmp_empty_mask);
            
            % Check the zy profile
            fig_hdl = figure;
            ax_1 = subplot(2,1,1);
            imshowpair(squeeze(ref_im(:, reg_vis_sec, :)).', squeeze(tmp_reg(:, reg_vis_sec, :)).');
            ax_1.Title.String = sprintf('Stack %d Displacement vector [%d, %d, %d]',...
                i, tmp_vec);
            ax_2 = subplot(2,1,2);
            imshowpair(squeeze(ref_im(:, reg_vis_sec, :)).', squeeze(tmp_reg_mv(:, reg_vis_sec, :)).');
            fprintf('Stack %d Displacement vector [%d, %d, %d]\n', i, tmp_vec);
            disp_vec(:, i, reg_ch) = tmp_vec;
            exp_im_reg{i, reg_ch} = tmp_reg_mv;
            exp_im_reg{i, non_reg_ch} = imtranslate(exp_im_reg{i, non_reg_ch}, tmp_vec);            
        end
    end
end

% im_diff = exp_im_avg{1} - exp_im_avg{2};
% im_size = size(im_diff);
%% Register the image before ablation to the image after first ablation 
% This is a bit confusing... 
reg_ch = 1;
non_reg_ch = 2;
exp_im_reg = exp_im_avg;
num_sets = size(exp_im_avg, 1);
disp_vec = zeros(3, num_sets, num_ch);
reg_vis_sec = 600;
if registration_Q
    mfft_th = 1e4;
    fft_search_range = [1, 1, 15];
    
    ref_im = exp_im_reg{2, reg_ch};
    tmp_reg = exp_im_reg{1, reg_ch};
    tmp_ref_mask = ref_im >= mfft_th;
    tmp_reg_mask = tmp_reg >= mfft_th;
    [tmp_vec, tmp2] = MaskedTranslationRegistration(ref_im, tmp_reg, ...
        tmp_ref_mask, tmp_reg_mask, fft_search_range);
    fprintf('Displacement vector [%d, %d, %d], correlation %.2f\n', tmp_vec, tmp2);
    tmp_reg_mv = imtranslate(tmp_reg, tmp_vec);
%     tmp_empty_mask = (tmp_reg_mv == 0);
%     tmp_reg_mv(tmp_empty_mask) = tmp_reg(tmp_empty_mask);
    exp_im_reg{1, reg_ch} = tmp_reg_mv;
    exp_im_reg{1, non_reg_ch} = imtranslate(exp_im_reg{1, non_reg_ch}, tmp_vec);
    disp_vec(:, i, reg_ch) = tmp_vec;
end
%% Merge channels 
merge_reg = cell(num_sets, 1);
for i = 1 : num_sets
    merge_reg{i} = fun_merge_image_stacks(exp_im_reg(i, :), 'method', 'GraySkull', 'stretch_contrast_Q', true);
end
%% Check registration result
% Without registration
abl_idx = 1;
vis_ch = 1;
abl_z_pxl = 50;
vis_im_stack = fun_merge_image_stacks(exp_im_avg([0, 1] + abl_idx, vis_ch), 'method', 'rgb', 'stretch_contrast_Q', true);
vis_im_stack_zyx = permute(vis_im_stack, [4, 1, 3, 2]);
vis_im_stack_zxy = permute(vis_im_stack, [4, 2, 3, 1]);

vis_im_stack_r = fun_merge_image_stacks(exp_im_reg([0, 1] + abl_idx, vis_ch), 'method', 'rgb', 'stretch_contrast_Q', true);
vis_im_stack_r_zyx = permute(vis_im_stack_r, [4, 1, 3, 2]);
vis_im_stack_r_zxy = permute(vis_im_stack_r, [4, 2, 3, 1]);
    %% Check registration result
vis_sec = 314;
fig_hdl = figure;
ax_hdl = subplot(2,1,1);
image(ax_hdl, squeeze(vis_im_stack_zxy(:, :, :, vis_sec)));
ax_hdl.DataAspectRatio = [1,1,1];
l_hdl = line(ax_hdl, ax_hdl.XLim, [abl_z_pxl, abl_z_pxl], 'Color', 'w', 'LineStyle', '-.');
ax_2 = subplot(2,1,2);
image(ax_2, squeeze(vis_im_stack_r_zxy(:, :, :, vis_sec)));
ax_2.DataAspectRatio = [1,1,1];
ax_2.Title.String = sprintf('Translation [%d, %d, %d] \\mum', disp_vec(:, abl_idx + 1, 1));
l_hdl_2 = line(ax_2, ax_2.XLim, [abl_z_pxl, abl_z_pxl], 'Color', 'w', 'LineStyle', '-.');
    %% MIP - ZYX
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 1.5;
ax_hdl = subplot(2,1,1);
image(ax_hdl, squeeze(max(vis_im_stack_zyx, [], 4)));
ax_hdl.DataAspectRatio = [1,1,1];
l_hdl = line(ax_hdl, ax_hdl.XLim, [abl_z_pxl, abl_z_pxl], 'Color', 'w', 'LineStyle', '-.');
ax_2 = subplot(2,1,2);
image(ax_2, squeeze(max(vis_im_stack_r_zyx, [], 4)));
ax_2.DataAspectRatio = [1,1,1];
ax_2.Title.String = sprintf('Translation [%d, %d, %d] \\mum', disp_vec(:, 2, 1));
l_hdl_2 = line(ax_2, ax_2.XLim, [abl_z_pxl, abl_z_pxl], 'Color', 'w', 'LineStyle', '-.');
ax_hdl.Title.String = sprintf('700 \\mum MIP after %d \\mum ablation', abl_z_pxl);
ax_hdl.YLabel.String = sprintf('Z (\\mum)');
ax_2.XLabel.String = 'Y (\mum)';
% fig_fp = fullfile(exp_fig_folder, sprintf('%s_50um_abl_%d_zyx_700um_mip_ch%d.png', exp_fig_prefix, abl_idx, vis_ch));
% fun_print_image_in_several_formats(fig_hdl, fig_fp);
    %% MIP - ZXY
mip_range = 1 : 700;
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 1.5;
ax_hdl = subplot(2,1,1);
image(ax_hdl, squeeze(max(vis_im_stack_zxy(:, :, :, mip_range), [], 4)));
ax_hdl.DataAspectRatio = [1,1,1];
l_hdl = line(ax_hdl, ax_hdl.XLim, [abl_z_pxl, abl_z_pxl], 'Color', 'w', 'LineStyle', '-.');
ax_2 = subplot(2,1,2);
image(ax_2, squeeze(max(vis_im_stack_r_zxy(:, :, :, mip_range), [], 4)));
ax_2.DataAspectRatio = [1,1,1];
ax_2.Title.String = sprintf('Translation [%d, %d, %d] \\mum', disp_vec(:, 2, 1));
l_hdl_2 = line(ax_2, ax_2.XLim, [abl_z_pxl, abl_z_pxl], 'Color', 'w', 'LineStyle', '-.');
ax_hdl.Title.String = sprintf('700 \\mum MIP after %d \\mum ablation', abl_z_pxl);
ax_hdl.YLabel.String = sprintf('Z (\\mum)');
ax_2.XLabel.String = 'X (\mum)';
fig_fp = fullfile(exp_fig_folder, sprintf('%s_50um_abl_%d_zxy_700um_mip_ch%d.png', exp_fig_prefix, abl_idx, vis_ch));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
    %% Compare before and after registration 
vis_sec = 50;
fig_hdl = figure;
ax_1 = subplot(1,2,1);
imshowpair(exp_im_avg{1, 1}(:, :, vis_sec), exp_im_avg{2, 1}(:, :, vis_sec));
ax_1.Title.String = sprintf('Raw @ %d \\mum', vis_sec);
ax_2 = subplot(1,2,2);
imshowpair(exp_im_reg{1, 1}(:, :, vis_sec), exp_im_reg{2, 1}(:, :, vis_sec));
ax_2.Title.String = sprintf('After registration @ %d \\mum', vis_sec);
%% Side view, single stack
vis_im_zxy = permute(merge_reg{3}, [4, 2, 3, 1]);
vis_sec_list = 200 : 100 : 600;
num_vis = numel(vis_sec_list);

fig_hdl = figure;
fig_hdl.Position(4) = fig_hdl.Position(4) * 2;
% ax_hdl = axes(fig_hdl);
for i = 1 : num_vis
    vis_sec = vis_sec_list(i);
    ax_hdl = subplot(num_vis, 1, i);
    image(ax_hdl, vis_im_zxy(:, :, :, vis_sec));
    ax_hdl.DataAspectRatio = [1,1,1];
    hold(ax_hdl, 'on');
    rectangle(ax_hdl, 'Position', [100, 50, 500, 50], 'LineStyle', '--', 'EdgeColor', 'r');
    ax_hdl.XLabel.String = 'Z (\mum)';
    ax_hdl.YLabel.String = 'X (\mum)';
    ax_hdl.Title.String = sprintf('Y = %d \\mum', vis_sec);
end
% tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_abl_%d_aft_vis_zxy_secs.png',...
%     exp_fig_prefix, abl_idx));
% fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
%%
vis_ch = 2;
tmp_bbox_ctr = round(im_size(1:2) / 2);
tmp_bbox_ll_um = [500, 500];
tmp_bbox_ll = tmp_bbox_ll_um / pixel_size;
tmp_bbox_mmll = [tmp_bbox_ctr - tmp_bbox_ll/2, tmp_bbox_ll];

vis_sec_idx = [10, 30, 55, 60, 65];
num_sec = numel(vis_sec_idx);
fig_hdl = figure;
fig_hdl.Position(3) = fig_hdl.Position(3) * 4;
for i = 1 : num_sec
    ax_1 = subplot(2, num_sec, sub2ind([num_sec, 2], i, 1));
    ax_2 = subplot(2, num_sec, sub2ind([num_sec, 2], i, 2));
    imagesc(ax_1, exp_im_reg{1, vis_ch}(:, :, vis_sec_idx(i)));
    imagesc(ax_2, exp_im_reg{2, vis_ch}(:, :, vis_sec_idx(i)));
    colormap(ax_1, 'gray');
    colormap(ax_2, 'gray');
    ax_1.CLim = [0, 60000];
    ax_2.CLim = [0, 60000];
    [ax_1.DataAspectRatio, ax_2.DataAspectRatio] = deal([1,1,1]);
    if i <= 3
        hold(ax_2, 'on');
        rectangle(ax_2, 'Position', tmp_bbox_mmll, 'LineStyle', '--', 'EdgeColor', 'r');
    end
    if i == 1
        ax_1.YLabel.String = 'Before';
        ax_2.YLabel.String = 'After';
        ax_1.XTick = [];
        [ax_1.YTick, ax_2.YTick] = deal([]);
        ax_2.XTickLabel = arrayfun(@(x) num2str(x * pixel_size), ax_2.XTick, 'UniformOutput', false);
        ax_2.XLabel.String = 'X (\mum)';
    else
        [ax_1.XTick, ax_2.XTick] = deal([]);
        [ax_1.YTick, ax_2.YTick] = deal([]);
    end
    ax_1.Title.String = sprintf('%d \\mum', vis_sec_idx(i));
end

exp_fig_folder = fullfile(vis_root_folder, site_folder);
exp_fig_prefix = sprintf('%s_%s', vis_file_prefix, site_folder);
tmp_fig_fp = fullfile(exp_fig_folder, sprintf('%s_yx_abl_1_vs_2_ch%d.png',...
    exp_fig_prefix, vis_ch));
fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
%%

