function im = record_construction_state()
% Take a picture and save to folder
user_folder_path = fullfile(char(java.lang.System.getProperty('user.home')), ...
    'Pictures', 'WBIM');
file_name = sprintf('WBIM_construction_record_%s.jpg', datestr(datetime('now'), 'yyyymmdd_HHMMSS'));
fn = fullfile(user_folder_path, file_name);
cam = webcam;
im = snapshot(cam);

if ~isfolder(user_folder_path)
    mkdir(user_folder_path);
end
imwrite(im, fn);
end