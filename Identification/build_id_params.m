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
params.late_time = 80.0;

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

end