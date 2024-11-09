DataManager = WBIMFileManager;
data_folder = '/net/birdstore/Vessel/LifeCanvas/Moreno_DK10523-001_15x_Fluorescein-Autofluorescence/488nm_Fluorescein';
[parent_folder, data_folder_name] = fileparts(data_folder);
im_list = dir(fullfile(data_folder, '*.tif'));
voxel_size_um = [0.41, 0.41, 1];
num_im = numel(im_list);
% Resize to 1 um isotropic resolution 
%%
file_name_label = arrayfun(@(x) strsplit(x.name(1:end-4), '_'), im_list, 'UniformOutput', false);
file_name_label = cat(1, file_name_label{:});
file_name_label = file_name_label(:, 3);
file_name_label = cellfun(@(x) str2num(x), file_name_label);
assert(issorted(file_name_label, 'ascend'), 'Filenames are not sorted');
%%
im_size = [9230 7422];
im_size_um = round(im_size .* voxel_size_um(1:2));
stack_size = [im_size_um, num_im];
im_stack = zeros(stack_size, 'uint16');
t_tic = tic;
for i = 1 : num_im
    tmp_im = DataManager.load_single_tiff(fullfile(im_list(i).folder, ...
        im_list(i).name));
    im_stack(:, :, i) = imresize(tmp_im, im_size_um);
    if mod(i, 100) == 0
       fprintf("Finish loading %d (%.3f) images. Elapsed time is %.2f second\n",...
           i, i/num_im, toc(t_tic));
    end
end
%%
im_stack_zyx = permute(im_stack, [3, 1, 2]);
im_stack_zyx = fun_stretch_contrast(im_stack_zyx, 0.005, 0.99);
implay(im_stack_zyx);
fig_fp = fullfile(parent_folder, sprintf('%s_sagittal.avi', data_folder_name));
fun_vis_write_stack_to_avi(im2uint8(im_stack_zyx), fig_fp);
%%
im_stack_zxy = permute(im_stack_zyx, [1, 3, 2]);
implay(im_stack_zxy);
fig_fp = fullfile(parent_folder, sprintf('%s_coronal.avi', data_folder_name));
fun_vis_write_stack_to_avi(im2uint8(im_stack_zxy), fig_fp);
