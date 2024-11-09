clc;clear;close all;
DataManager = WBIMFileManager;
exp_group = 'TestSample';
exp_name = 'WBIM20230414002_20230508';
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization', ...
    'abl_roi_detection');
im_root_folder = '/net/birdstore/Vessel/WBIM/Acquisition/TestSample/WBIM20230414002_20230508/00001/Explore';
%%
tissue_fps = dir(fullfile(im_root_folder, 'Tissue*'));
num_tissue_im = numel(tissue_fps);
tissue_im_cell = cell(num_tissue_im, 1);
num_sec = 200;
for i = 1 : num_tissue_im
    tmp_im = DataManager.load_single_tiff(fullfile(tissue_fps(i).folder, ...
        tissue_fps(i).name));
    tmp_num_ch = size(tmp_im, 3) / num_sec;
    tmp_ch_im = arrayfun(@(x) tmp_im(:, :, x:tmp_num_ch:end), 1:tmp_num_ch, ...
        'UniformOutput', false);
    % Median filter to smooth the image
    tmp_ch_im = cellfun(@medfilt3, tmp_ch_im, 'UniformOutput', false);
    tissue_im_cell{i} = tmp_ch_im;    
end
%% Chromatic aberration
% Merge image from different wavelength - same channel 
vsl_ch_im = cellfun(@(x) x{1},  tissue_im_cell(4:6), 'UniformOutput', false);
vsl_ch_im_mip = cellfun(@(x) max(x, [], 3), vsl_ch_im, 'UniformOutput', false);
vsl_ch_im_rgb = cat(4, vsl_ch_im{:});
vsl_ch_im_rgb = permute(vsl_ch_im_rgb, [1,2,4,3]);
implay(uint16(vsl_ch_im_rgb)*2)
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = nexttile();
imagesc(vsl_ch_im_mip{1} - vsl_ch_im_mip{2});
ax_hdl.DataAspectRatio = [1,1,1];
c_hdl = colorbar(ax_hdl);
ax_hdl.Title.String = 'I_{1030} - I_{1060}';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
colormap(ax_hdl, 'jet');
fig_fp = fullfile(vis_folder, sprintf('Tissue_1030nm_-_1060nm.png'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

% imshowpair(vsl_ch_im_mip{2}, vsl_ch_im_mip{3});
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = nexttile();
imagesc(rescale(abs(vsl_ch_im_mip{3} - vsl_ch_im_mip{1})) - rescale(abs(vsl_ch_im_mip{1} - vsl_ch_im_mip{2})));
ax_hdl.DataAspectRatio = [1,1,1];
c_hdl = colorbar(ax_hdl);
ax_hdl.Title.String = '|I_{980} - I_{1030}| - |I_{1030} - I_{1060}|';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
colormap(ax_hdl, 'jet');
%% Tissue & tdTomato 
ts_ch_im = cellfun(@(x) x{2},  tissue_im_cell(4:6), 'UniformOutput', false);
ts_ch_im_mip = cellfun(@(x) max(x, [], 3), ts_ch_im, 'UniformOutput', false);
im_name = {tissue_fps(4:6).name};
im_name = cellfun(@(x) strsplit(x, '_'), im_name, 'UniformOutput', false);
im_name = cellfun(@(x) x{2}, im_name, 'UniformOutput', false);
fig_hdl = figure;
for i = [3,1,2]
    ax_hdl = nexttile();
    imagesc(ax_hdl, ts_ch_im_mip{i});
    ax_hdl.DataAspectRatio = [1,1,1];
    colorbar;
    ax_hdl.Title.String = im_name{i};
end