data_folder = 'H:\WBIM\Acquisition\PSF\20231222';
data_file = 'ObjNA095_0025pwr_050hwp_35mag_10fps_025umz_900nm_00001.tif';
data = DataManager.load_single_tiff(fullfile(data_folder, data_file));

%%
num_avg = 10;
data_avg = zeros([size(data, [1,2]), size(data, 3) / num_avg], 'like', ...
    data);
for i = 1 : size(data_avg, 3)
    data_avg(:, :, i) = mean(medfilt3(data(:, :, (i-1) * num_avg + 1 : ...
        i * num_avg), [3,3,1]), 3);
end
data_avg = uint16(data_avg) * 2;
%%
ctr_xy = [317, 175];
z_prof = data_avg(ctr_xy(2), ctr_xy(1), :);
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
plot(squeeze(z_prof))
% 1, 2 were using Pantong's objective
% 3 was using 8 mm 0.95 na objective