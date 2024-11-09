DataManager = WBIMFileManager;
exp_group = 'Ablation_test';
exp_name = 'WBIMGel_20220901';
data_folder = DataManager.fp_layer(exp_group, exp_name, 1);
image_folder = fullfile(data_folder, 'Explore');
vis_folder = fullfile(DataManager.fp_experiment(exp_group, exp_name), 'visualization');
im_files = dir(fullfile(image_folder, '*.tif'));
%%
im_file_names = {im_files.name}';
im_file_names_d = cellfun(@(x) strsplit(x(1:end-4), '_'), im_file_names, 'UniformOutput', false);
im_file_names_d = cat(1, im_file_names_d{:});
im_site_index = cellfun(@(x) str2double(x), im_file_names_d(:, 2));
im_site_im_idx = cellfun(@(x) str2double(x), im_file_names_d(:, 3));
[bin_idx, bin_val] = fun_bin_data_to_idx_list(im_site_index);
file_info = struct;
file_info.site_index = bin_val;
file_info.filename = cellfun(@(x) im_file_names(x), bin_idx, 'UniformOutput', false);
file_info.filepath = cellfun(@(x) fullfile(image_folder, im_file_names(x)), bin_idx, 'UniformOutput', false);
file_info.time = arrayfun(@(x) datestr(x.datenum, 'YYYYmmDD_hhMMss'), im_files, 'UniformOutput', false);
%% Abl_settings
abl_info_file = dir(fullfile(data_folder, 'ablation', '*.mat'));
abl_info = arrayfun(@(x) load(fullfile(x.folder, x.name)), abl_info_file, 'UniformOutput', false);
abl_info = cat(1, abl_info{:});
%%
site_index = 3;
site_file = cellfun(@(x) DataManager.load_single_tiff(x), file_info.filepath{site_index}, 'UniformOutput', false);
site_file = cellfun(@(x) uint16(x) * 2, site_file, 'UniformOutput', false);

site_file_sm = cellfun(@(x) medfilt3(x), site_file, 'UniformOutput', false);
volumeViewer(site_file_sm{1} - site_file_sm{2})

