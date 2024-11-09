file_path = mfilename('fullpath');
[file_folder, ~, ~] = fileparts(file_path);
if isunix
    path_sep = ':';
elseif ispc
    path_sep = ';';
end
folders_path = strsplit(genpath(file_folder), path_sep);
is_git_Q = contains(folders_path, '.git');
folders_path = strjoin(folders_path(~is_git_Q), ';');
addpath(folders_path);
DataManager = WBIMFileManager;
% addpath(genpath(DataManager.SCRIPT_PATH));
[github_folder, ~, ~] = fileparts(DataManager.SCRIPT_PATH);
util_folder = fullfile(github_folder, 'MUtil');
switch DataManager.HOSTNAME
    case 'ThinkPadP16'
        scan_image_fp = 'C:\Program Files\Vidrio\SI2021.0.0_2024-02-08-130623_3baa4b74be';
        addpath(genpath(scan_image_fp));
%         fiji_fp = 'C:\Users\xij072\Software\Fiji.app\scripts';
%         addpath(fiji_fp);
    case 'pia'
        fiji_fp = 'C:\Users\dklab\Documents\GitHub\Fiji\scripts';
        addpath(fiji_fp);
        % Might need to change this: C:\Program Files\Vidrio\Old_version\SI2021.0.0_2021-07-30-144329_3baa4b74be\+scanimage\+components\+scan2d\+rggscan\Control.m
        % Add Digital refractometer
%         RF_fp = 'C:\Users\dklab\Documents\GitHub\Digital_Refractometer';
%         addpath(RF_fp);
    case {'bird', 'bird.dk.ucsd.edu'}
        lens_def_fp = fullfile(github_folder, 'Lens_Deformation');
        addpath(genpath(lens_def_fp));             
    otherwise 
        fprintf('Unrecognized machine. Unable to add ScanImage folder.\n');
end
if isfolder(util_folder)
    folders_path = strsplit(genpath(util_folder), path_sep);
    is_git_Q = contains(folders_path, '.git');
    folders_path = strjoin(folders_path(~is_git_Q), path_sep);
    addpath(folders_path);
end
clc;clear;close all;