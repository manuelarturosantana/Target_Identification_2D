function params = build_id_params(varargin)
%BUILD_ID_PARAMS  Build parameter struct for target identification.
%
%   params = BUILD_ID_PARAMS() returns a struct containing default
%   parameters.
%
%   params = BUILD_ID_PARAMS(name, value, ...) overrides default parameter
%   values using name-value pairs.
%
%   Only parameters that already exist as default fields may be overridden.
%   An error is raised if an unknown parameter name is supplied.
%
%   Example:
%       params = build_id_params("late_time", 100);
%
%   Current default fields:
%       late_time   Late-time cutoff used for classification.

% Default values
params = struct();

% paths
params.input_folder = 'local/data_two_cavity_constellations_2d';
params.output_path = ['Identification/classification_results/' ...
                      'constellation_results.mat'];

% classification parameters
params.late_time = 80.0;
params.noise_sigmas = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
params.svd_tol = 1e-12;
params.random_seed = 2026;

% Other options are averaged, single_objected_averaged
params.snr_report_type = "representative_object"; 
% averaged:              means average_snr overall obstacles on a per_angle basis.
%                        variance is also reported in this case
% single_object_averaged: Average over one object. Variance is not
%                         reported. In this case the field snr_object_name is used to 
%                         choose which object should be averaged.

% Representative signal used for SNR reporting.
% Defaults mean: first object/configuration, first receiver.
params.snr_object_name = "";
params.snr_receiver_angle = [];



% additional parameters
params.show_progress = false;

% Check name-value pair structure
if mod(numel(varargin), 2) ~= 0
    error("build_id_params:InvalidInput", ...
          "Inputs must be name-value pairs.");
end

% Override defaults
for k = 1:2:numel(varargin)
    name = varargin{k};
    value = varargin{k + 1};

    if ~ischar(name) && ~isstring(name)
        error("build_id_params:InvalidName", ...
              "Parameter names must be strings or character arrays.");
    end

    name = char(name);

    if ~isfield(params, name)
        error("build_id_params:UnknownParameter", ...
              "Unknown parameter '%s'.", name);
    end

    params.(name) = value;
end

% Validate snr_report_type

stype = char(params.snr_report_type);
% canonical options
valid = {'averaged', 'single_object_averaged', ...
         'representative_object'};
if ~any(strcmp(stype, valid))
    error("build_id_params:UnknownParameterValue", ...
          "Unknown snr_report_type '%s'.", ...
          stype);
end
end