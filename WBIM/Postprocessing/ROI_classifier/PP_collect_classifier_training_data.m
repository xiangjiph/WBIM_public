DataManager = WBIMFileManager;
exp_group_name = 'TestSample';
exp_name = 'WBIM20220610';
exp_root = fullfile(DataManager.fp_acquisition_disk, ...
    DataManager.fpr_experiment(exp_group_name, exp_name));
%% Get all the explored MIP
mip_dir = dir(fullfile(exp_root, '**', 'Explore', '**/mip_CH1.tif'));
num_mip = numel(mip_dir);
mip_data = cell(1, num_mip);
for i = 1 : num_mip
    tmp_fp = fullfile(mip_dir(i).folder, mip_dir(i).name);
    mip_data{i} = DataManager.load_single_tiff(tmp_fp);
end
%%
im_labler = ImageLabelingApp(mip_data);
im_label = im_labler.data_label;
%% Split images into 100 x 100 um2 patches





%% Compute features
im_stat = cellfun(@WBIMAcqPostProcessing.analyze_single_tile_mip, mip_data);
is_fg_or_bg_Q = (im_label == 0 | im_label == 1);
bg_vs_fg_data = struct;
bg_vs_fg_data.features = struct2table(im_stat);
bg_vs_fg_data.label = logical(im_label(is_fg_or_bg_Q)).';
selected_field_name = {'mean_mft', 'std_mft', 'max_mft', 'mean_rf', 'std_rf', 'max_rf', ...
    'frac_above_half_max_rf', 'frac_saturated_rf', 'nmpdf', 'mpdf'};

bg_vs_fg_classifier = fun_learning_get_classifier(bg_vs_fg_data, selected_field_name);

%% Save data
training_data = struct;
training_data.source_file = mip_dir;
% training_data.data = mip_data;
training_data.label = im_label;
training_data.label_name = {'Background', 'Foreground', 'Boundary'};
training_data.exp_group = exp_group_name;
training_data.experiment_name = exp_name;
training_data.data_description = 'Explore mode MIP CH1';
training_data.filepath = fullfile(exp_root, sprintf('%s_%s_classifier_training_data_%s.mat', ...
    exp_group_name, exp_name, datestr(now, 'yyyymmdd')));
save(training_data.filepath, '-struct', 'training_data');