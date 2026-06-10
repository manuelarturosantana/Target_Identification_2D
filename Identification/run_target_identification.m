% run_identification.m
%
% Run GLRT target identification from precomputed pole/time-domain data.

clear
clc

% ---- user-defined parameters ----
% paths
input_folder = 'local/data_two_cavity_constellations_2d';

output_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'constellation_results.mat');

% classification params
late_time = 80.0;
sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];

% SNR based on representative object
% If snr_object_name is empty, the first object is used.
% If snr_receiver_angle is empty, the first receiver is used.
snr_object_name = 'single_circ_o050';
snr_receiver_angle = 3*pi/2;

% build params struct
params = build_id_params( ...
    'input_folder', input_folder, ...
    'output_path', output_path, ...
    'late_time', late_time, ...
    'noise_sigmas', sigmas, ...
    'snr_object_name', snr_object_name, ...
    'snr_receiver_angle', snr_receiver_angle, ...
    'show_progress', true);

% load objects and spectral data (poles, time-domain signals)
data = load_pole_data(params);

% cache matrices and SVDs used in GLRT computation (speeds up script)
model = build_glrt_model(data, params);

% for each time-domain signal, compute GLRT statistic on library, classify
results = classify_glrt_data(data, model, params);

% save classification results
output_folder = fileparts(params.output_path);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

save(params.output_path, 'results');

fprintf('Saved classification results to:\n  %s\n', params.output_path);