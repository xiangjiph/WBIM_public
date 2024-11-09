exp_root_folder = 'G:\WBIM\Acquisition\TestSample\WBIM20230124001_20230301\00001\ablation\Ablation_test';
exp_folder = fullfile(exp_root_folder, 'Site4');
exp_files = dir(fullfile(exp_folder, '*.tif'));
%%
num_ch = 2;
num_files = numel(exp_files);
exp_data = cell(num_files, 1);
for i = 1 : num_files
    tmp_data_fp = fullfile(exp_files(i).folder, exp_files(i).name);
    tmp_data = uint16(DataManager.load_single_tiff(tmp_data_fp))*2;
    tmp_stack_size = size(tmp_data);
    tmp_stack_size(3) = tmp_stack_size(3) / num_ch;
    tmp_data_rgb = zeros([tmp_stack_size, 3], 'uint8');
    for j = 1 : num_ch
        tmp_data_rgb(:, :, :, j) = im2uint8(medfilt3(tmp_data(:, :, j:num_ch:end)));
    end
    exp_data{i} = permute(tmp_data_rgb, [1,2,4,3]);
end
%%
implay(exp_data{1})
implay(exp_data{3})

implay(permute(exp_data{1}, [4, 1, 3, 2]));