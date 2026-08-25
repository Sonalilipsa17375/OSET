% An example script for extracting features from ECG signals
% Authors: 
% Seyedeh Somayyeh Mousavi
% Reza Sameni
% bmemousavi@gmail.com
% Aug 2026
% Emory University, Georgia, USA
%% ====================================================================
clc
clear
close all
%% ====================================================================
% Add OSET package to path
path_function_1 = '../oset_ecg_feature_extraction/';
path_function_2 = '../../../../OSET/';
addpath(genpath(path_function_1));
addpath(genpath(path_function_2));
%% ====================================================================
% Input and output paths
input_path = 'sample-ecg';
output_path = 'features_csv_files';
extract_ecg_features(input_path, output_path)
disp('Finish')
%% ====================================================================
% Remove path
rmpath(genpath(path_function_1));
rmpath(genpath(path_function_2));