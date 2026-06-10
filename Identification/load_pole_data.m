function data = load_pole_data(params)
%LOAD_POLE_DATA  Load pole and time-domain data for target identification.
%
%   data = LOAD_POLE_DATA(params) loads all .mat files in
%   params.input_folder. The returned struct array has one entry per target
%   configuration.
%
%   Fields:
%       data(k).name             configuration name
%       data(k).file_name        .mat file name
%       data(k).file_path        full path to original .mat file
%       data(k).poles            complex poles
%       data(k).ts               time grid
%       data(k).signals          receiver time-domain signals
%       data(k).receiver_angles  receiver angles, if available

    input_folder = params.input_folder;

    files = dir(fullfile(input_folder, '*.mat'));

    if isempty(files)
        error('load_pole_data:NoFiles', ...
              'No .mat files found in folder: %s', input_folder);
    end

    % Sort by file name so object indices are reproducible.
    [~, idx] = sort({files.name});
    files = files(idx);

    data = struct([]);

    for k = 1:numel(files)
        file_path = fullfile(files(k).folder, files(k).name);

        S = load(file_path, ...
            'config_name', ...
            'pols', ...
            'ts', ...
            'uff_all', ...
            'receiver_angles');

        % Use config_name if present. Otherwise infer name from file name.
        if isfield(S, 'config_name')
            name = string(S.config_name);
        else
            [~, stem, ~] = fileparts(files(k).name);
            stem = string(stem);

            tokens = regexp(stem, '^\d+_(.*)$', 'tokens', 'once');

            if isempty(tokens)
                name = stem;
            else
                name = string(tokens{1});
            end
        end

        % Required fields.
        if ~isfield(S, 'pols')
            error('load_pole_data:MissingPoles', ...
                  'File %s is missing variable "pols".', file_path);
        end

        if ~isfield(S, 'ts')
            error('load_pole_data:MissingTimeGrid', ...
                  'File %s is missing variable "ts".', file_path);
        end

        if ~isfield(S, 'uff_all')
            error('load_pole_data:MissingSignals', ...
                  'File %s is missing variable "uff_all".', file_path);
        end

        signals = S.uff_all;

        % Convention: one row per receiver, one column per time sample.
        if isvector(signals)
            signals = signals(:).';
        end

        data(k).name = name;
        data(k).file_name = string(files(k).name);
        data(k).file_path = string(file_path);
        data(k).poles = S.pols(:);
        data(k).ts = S.ts(:);
        data(k).signals = signals;

        if isfield(S, 'receiver_angles')
            data(k).receiver_angles = S.receiver_angles(:);
        else
            data(k).receiver_angles = [];
        end
    end

    % Basic consistency checks.
    ts0 = data(1).ts;
    nt = numel(ts0);

    for k = 1:numel(data)
        if numel(data(k).ts) ~= nt || ...
                max(abs(data(k).ts - ts0)) > 1e-12
            error('load_pole_data:TimeGridMismatch', ...
                  'Time grid mismatch in file: %s', data(k).file_path);
        end

        if size(data(k).signals, 2) ~= nt
            error('load_pole_data:SignalLengthMismatch', ...
                  ['Signal length does not match time grid in file:\n' ...
                   '  %s\n' ...
                   'Expected %d time samples, got %d.'], ...
                  data(k).file_path, nt, size(data(k).signals, 2));
        end

        if ~isempty(data(k).receiver_angles) && ...
                numel(data(k).receiver_angles) ~= size(data(k).signals, 1)
            error('load_pole_data:ReceiverMismatch', ...
                  ['Number of receiver angles does not match number ' ...
                   'of signals in file:\n  %s'], ...
                  data(k).file_path);
        end
    end

end