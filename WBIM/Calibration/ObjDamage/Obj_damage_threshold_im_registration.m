fp_im_1 = 'D:\My Drive\Microscope\Images\Test objective\20220214\Obj_20220214_before_8x_2.tif';
fp_im_2 = 'D:\My Drive\Microscope\Images\Test objective\20220214\Obj_20220214_aft_100_1hr_8x_3.tif';

im_1 = rescale(imread(fp_im_1));
im_2 = rescale(imread(fp_im_2));
%%
[im_1_g, ~] = imgradient(im_1);
im_1_g_sum = sum(im_1_g, 1);

[im_2_g, ~] = imgradient(im_2);
im_2_g_sum = sum(im_2_g, 1);

plot(im_1_g_sum)
hold on
plot(im_2_g_sum)
%% Crop the image
im_1_c = im_1(:, 1050 : 2150);
im_2_c = im_2(:, 750 : 1900);
[im_1_c_g, ~] = imgradient(im_1_c);
[im_2_c_g, ~] = imgradient(im_2_c);
%% 
[optimizer, metric] = imregconfig('multimodal');
% optimizer.MaximumStepLength = 0.03;
reg_tform = imregtform(im_1_c_g, im_2_c_g, 'affine', optimizer, metric, ...
    'DisplayOptimization', true);
im_1_c_reg = imwarp(im_1_c, reg_tform, 'OutputView', imref2d(size(im_2_c)));
%%
imshowpair(im_2_c, im_1_c_reg);
%%
fig_hdl = figure;
subplot(1,2,1);
imagesc(im_1_c);
set(gca, 'DataAspectRatio', [1,1,1]);
subplot(1,2,2);
imagesc(im_2_c);
set(gca, 'DataAspectRatio', [1,1,1]);