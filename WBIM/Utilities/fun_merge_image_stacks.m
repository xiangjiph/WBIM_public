function im = fun_merge_image_stacks(im_cell, opts)

arguments
    im_cell (1, :) cell
    opts.method = 'rgb';
    opts.stretch_contrast_Q (1,1) logical = false;
    opts.output_dtype = 'uint8';
end
num_ch = numel(im_cell);

if opts.stretch_contrast_Q
    im_cell = cellfun(@fun_stretch_contrast, im_cell, 'UniformOutput', false);
end
% Convert datatype
if ~strcmp(opts.output_dtype, im_cell{1})
    switch opts.output_dtype
        case 'uint8'
            cast_fun = @im2uint8;
        case 'uint16'
            cast_fun = @im2uint16;
        otherwise
            error('To be implemented');
    end
    im_cell = cellfun(cast_fun, im_cell, 'UniformOutput', false);
end
im_size = size(im_cell{1});
if numel(im_size) == 2
    im_size = [im_size, 1];
end
im = zeros([im_size, 3], opts.output_dtype);
switch opts.method
    case 'rgb'
        for i = 1 : num_ch
            im(:, :, :, i) = im_cell{i};
        end
    case 'VesselBone'
        im = repelem(im_cell{2}, 1, 1, 1, 3) * 0.6;
        im(:, :, :, 1) = max(im(:, :, :, 1), im_cell{1});      
    case 'GraySkull'
        % Assume the seccond channel is SHG
        im = repmat(im_cell{WBIMChannelName.SHG} .* 0.4, 1, 1, 1, 3);
        % Maximum 4 channels
        non_skull_ch = 1 : num_ch;
        non_skull_ch = non_skull_ch(non_skull_ch ~= WBIMChannelName.SHG);
        for i_ch = 1 : numel(non_skull_ch)
            im(:, :, :, i_ch) = max(im(:, :, :, i_ch), ...
                im_cell{non_skull_ch(i_ch)});
        end
    otherwise
        error('To be implemented');
end
im = permute(im, [1,2,4,3]);
im = squeeze(im);
end