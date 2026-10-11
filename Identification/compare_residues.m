function diagnostics = compare_residues()
%COMPARE_RESIDUES  Compare stored and time-domain fitted residues.
%
% Residues are fitted with the Data_Generation pole convention
%
%                 q(t) = sum_k r_k exp(-1i phi_k t).
%
% The complex time-domain data are retained throughout the fit.
% The script compares raw residues, relative modal powers, complex values,
% late-time reconstructions, and angle-dependent fitting errors.

% User parameters
input_folder = 'local/plane_data_360';
object_name = "twinjet_plane1";
receiver_angle = 3 * pi / 2;
late_time = 80.0;
fit_method = "least_squares";       % "ridge" or "least_squares"
ridge = 1e-10;
noise_sigma = 0;
random_seed = 1;
output_folder = fullfile('Identification', 'residue_diagnostics');

raw = load_selected_object(input_folder, object_name);
validate_object(raw);

ts = raw.ts(:);
poles = raw.pols(:);
angles = raw.receiver_angles(:);
raw_spatial_residues = raw.r_res_all;
modal_factor = waveform_residue_factor(raw, poles);
original_residues = raw_spatial_residues .* modal_factor.';
late_idx = find(ts >= late_time, 1, 'first');

if isempty(late_idx)
    error('compare_residues:InvalidLateTime', ...
        'late_time exceeds the final time sample.');
end

tt = ts(late_idx:end);
Z = pole_basis(tt, poles);
clean_signals = raw.uff_all(:, late_idx:end);

rng(random_seed);
noise = (noise_sigma / sqrt(2)) * (randn(size(clean_signals)) ...
    + 1i * randn(size(clean_signals)));
noisy_signals = clean_signals + noise;
coefficients = fit_coefficients(Z, noisy_signals.', fit_method, ridge);

fitted_residues = coefficients.';
fitted_signals = (Z * coefficients).';

original_signals = (Z * original_residues.').';

% A complex scale reveals whether the stored residues have the right modal
% pattern but use a different global amplitude and phase at a given angle.
scale_numerator = sum(conj(original_signals) .* clean_signals, 2);
scale_denominator = sum(abs(original_signals).^2, 2);
original_scale = scale_numerator ./ max(scale_denominator, eps);
scaled_original_residues = original_residues .* original_scale;
scaled_original_signals = original_signals .* original_scale;

original_power = normalize_rows(abs(original_residues).^2);
fitted_power = normalize_rows(abs(fitted_residues).^2);

residue_error = row_relative_error(fitted_residues, original_residues);
scaled_residue_error = row_relative_error( ...
    fitted_residues, scaled_original_residues);
power_hellinger = 0.5 * sum( ...
    (sqrt(original_power) - sqrt(fitted_power)).^2, 2);

original_signal_error = row_relative_error( ...
    original_signals, clean_signals);
scaled_original_signal_error = row_relative_error( ...
    scaled_original_signals, clean_signals);
fitted_signal_error = row_relative_error(fitted_signals, noisy_signals);

selected_idx = nearest_angle_index(angles, receiver_angle);
selected_angle = angles(selected_idx);
object_label = string(raw.config_name);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

object_tag = regexprep(char(object_label), '[^A-Za-z0-9_-]', '_');
object_tag = string(object_tag);
selected_figure_path = fullfile(output_folder, ...
    object_tag + "_selected_residue_comparison.pdf");
summary_figure_path = fullfile(output_folder, ...
    object_tag + "_residue_error_by_angle.pdf");
output_path = fullfile(output_folder, ...
    object_tag + "_residue_diagnostics.mat");

make_selected_figure(tt, poles, clean_signals, original_signals, ...
    scaled_original_signals, fitted_signals, original_residues, ...
    fitted_residues, original_power, fitted_power, selected_idx, ...
    selected_angle, object_label, selected_figure_path);

make_angle_summary(angles, residue_error, scaled_residue_error, ...
    power_hellinger, original_signal_error, ...
    scaled_original_signal_error, fitted_signal_error, object_label, ...
    summary_figure_path);

diagnostics = struct();
diagnostics.object_name = object_label;
diagnostics.file_path = string(raw.file_path);
diagnostics.fit_method = fit_method;
diagnostics.ridge = ridge;
diagnostics.noise_sigma = noise_sigma;
diagnostics.late_time = late_time;
diagnostics.poles = poles;
diagnostics.receiver_angles = angles;
diagnostics.raw_spatial_residues = raw_spatial_residues;
diagnostics.waveform_residue_factor = modal_factor;
diagnostics.original_residues = original_residues;
diagnostics.fitted_residues = fitted_residues;
diagnostics.scaled_original_residues = scaled_original_residues;
diagnostics.original_scale = original_scale;
diagnostics.original_relative_power = original_power;
diagnostics.fitted_relative_power = fitted_power;
diagnostics.residue_relative_error = residue_error;
diagnostics.scaled_residue_relative_error = scaled_residue_error;
diagnostics.power_hellinger_distance = power_hellinger;
diagnostics.original_signal_relative_error = original_signal_error;
diagnostics.scaled_original_signal_relative_error = ...
    scaled_original_signal_error;
diagnostics.fitted_signal_relative_error = fitted_signal_error;
diagnostics.basis_condition_number = cond(Z);
diagnostics.selected_receiver_idx = selected_idx;
diagnostics.selected_receiver_angle = selected_angle;
diagnostics.selected_clean_signal = clean_signals(selected_idx, :).';
diagnostics.selected_original_reconstruction = ...
    original_signals(selected_idx, :).';
diagnostics.selected_scaled_original_reconstruction = ...
    scaled_original_signals(selected_idx, :).';
diagnostics.selected_fitted_reconstruction = ...
    fitted_signals(selected_idx, :).';

save(output_path, 'diagnostics');

fprintf('Object: %s\n', object_label);
fprintf('Selected receiver angle: %.6f rad\n', selected_angle);
fprintf('Basis condition number: %.3e\n', diagnostics.basis_condition_number);
fprintf('Waveform residue relative error: %.3e\n', ...
    residue_error(selected_idx));
fprintf('Scale-aligned residue relative error: %.3e\n', ...
    scaled_residue_error(selected_idx));
fprintf('Relative-power Hellinger distance: %.3e\n', ...
    power_hellinger(selected_idx));
fprintf('Stored-residue signal error: %.3e\n', ...
    original_signal_error(selected_idx));
fprintf('Scale-aligned stored-residue signal error: %.3e\n', ...
    scaled_original_signal_error(selected_idx));
fprintf('Fitted signal error: %.3e\n', fitted_signal_error(selected_idx));
fprintf('Saved diagnostics to:\n  %s\n', output_path);
fprintf('Saved figures to:\n  %s\n  %s\n', ...
    selected_figure_path, summary_figure_path);
end


function raw = load_selected_object(input_folder, object_name)
%LOAD_SELECTED_OBJECT  Load a MAT file by config_name or filename stem.

files = dir(fullfile(input_folder, '*.mat'));

if isempty(files)
    error('load_selected_object:NoFiles', ...
        'No MAT files were found in %s.', input_folder);
end

target = string(object_name);

for i = 1:numel(files)
    file_path = fullfile(files(i).folder, files(i).name);
    metadata = load(file_path, 'config_name');
    [~, file_stem] = fileparts(files(i).name);

    name_matches = isfield(metadata, 'config_name') ...
        && strcmpi(string(metadata.config_name), target);
    file_matches = strcmpi(string(file_stem), target);

    if name_matches || file_matches
        raw = load(file_path);
        raw.file_path = file_path;
        return
    end
end

error('load_selected_object:ObjectNotFound', ...
    'No object named %s was found in %s.', target, input_folder);
end


function validate_object(raw)
%VALIDATE_OBJECT  Check fields and dimensions required by the diagnostic.

required = {'config_name', 'pols', 'r_res_all', ...
    'receiver_angles', 'ts', 'uff_all', 'psp'};

for i = 1:numel(required)
    if ~isfield(raw, required{i})
        error('validate_object:MissingField', ...
            'The selected MAT file has no %s field.', required{i});
    end
end

nangles = numel(raw.receiver_angles);
npoles = numel(raw.pols);

if ~isequal(size(raw.r_res_all), [nangles, npoles])
    error('validate_object:ResidueSizeMismatch', ...
        'r_res_all must have size nreceiver-by-npole.');
end

if size(raw.uff_all, 1) ~= nangles || size(raw.uff_all, 2) ~= numel(raw.ts)
    error('validate_object:SignalSizeMismatch', ...
        'uff_all must have size nreceiver-by-ntime.');
end
end


function factor = waveform_residue_factor(raw, poles)
%WAVEFORM_RESIDUE_FACTOR  Convert stored spatial residues to exp(-i*p*t) coefficients.
%
% Data_Generation stores r_res_all from comp_r_res. For a Gaussian pulse,
% sum_SEM applies
%
%   -1i * B(p) * r_res * exp(-1i*p*(t-t0)).
%
% Moving the time shift into the coefficient gives the factor below.

psp = raw.psp;
required = {'w_0', 'sigmas', 't_0'};

for i = 1:numel(required)
    if ~isfield(psp, required{i})
        error('waveform_residue_factor:MissingParameter', ...
            'The saved psp struct has no %s field.', required{i});
    end
end

poles = poles(:);
pulse_factor = exp(-(poles - psp.w_0).^2 / psp.sigmas);
factor = -1i * pulse_factor .* exp(1i * poles * psp.t_0);

if isfield(psp, 'has_zero_freq') && psp.has_zero_freq
    factor = 2 * factor;
end
end


function Z = pole_basis(times, poles)
%POLE_BASIS  Build exp(-1i * phi_k * t) for every pole phi_k.

Z = exp(-1i * times(:) * poles(:).');
end


function coefficients = fit_coefficients(Z, signals, fit_method, ridge)
%FIT_COEFFICIENTS  Fit modal coefficients for all receiver signals.

if strcmp(fit_method, "ridge")
    coefficients = (Z' * Z + ridge * eye(size(Z, 2))) \ (Z' * signals);
elseif strcmp(fit_method, "least_squares")
    coefficients = Z \ signals;
else
    error('fit_coefficients:UnknownMethod', ...
        'fit_method must be "ridge" or "least_squares".');
end
end


function values = normalize_rows(values)
%NORMALIZE_ROWS  Normalize nonnegative rows to sum to one.

values = values ./ max(sum(values, 2), eps);
end


function error_value = row_relative_error(approximation, reference)
%ROW_RELATIVE_ERROR  Compute one relative two-norm error per row.

error_value = vecnorm(approximation - reference, 2, 2) ...
    ./ max(vecnorm(reference, 2, 2), eps);
end


function idx = nearest_angle_index(angles, angle)
%NEAREST_ANGLE_INDEX  Find the nearest direction on the periodic circle.

distance = abs(mod(angles(:) - angle + pi, 2 * pi) - pi);
[~, idx] = min(distance);
end


function make_selected_figure(times, poles, clean_signals, original_signals, ...
        scaled_original_signals, fitted_signals, original_residues, ...
        fitted_residues, original_power, fitted_power, selected_idx, ...
        selected_angle, object_label, figure_path)
%MAKE_SELECTED_FIGURE  Plot detailed residue diagnostics at one angle.

npoles = numel(poles);
modes = (1:npoles).';
original = original_residues(selected_idx, :).';
fitted = fitted_residues(selected_idx, :).';

fig = figure('Color', 'w');
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile(layout)
bar(modes, [abs(original), abs(fitted)], 'grouped')
xlabel('Pole index')
ylabel('Residue magnitude')
legend('Stored', 'Fitted', 'Location', 'best')
grid on

nexttile(layout)
bar(modes, [original_power(selected_idx, :).', ...
    fitted_power(selected_idx, :).'], 'grouped')
xlabel('Pole index')
ylabel('Relative modal power')
legend('Stored', 'Fitted', 'Location', 'best')
ylim([0, 1])
grid on

nexttile(layout)
plot(real(original), imag(original), 'o', 'DisplayName', 'Stored')
hold on
plot(real(fitted), imag(fitted), 'x', 'DisplayName', 'Fitted')

for m = 1:npoles
    plot([real(original(m)), real(fitted(m))], ...
        [imag(original(m)), imag(fitted(m))], ':', ...
        'HandleVisibility', 'off')
    text(real(original(m)), imag(original(m)), sprintf(' %d', m))
end

xlabel('Real residue')
ylabel('Imaginary residue')
axis equal
grid on
legend('Location', 'best')

nexttile(layout)
plot(times, real(clean_signals(selected_idx, :)), ...
    'DisplayName', 'Time-domain signal')
hold on
plot(times, real(original_signals(selected_idx, :)), '--', ...
    'DisplayName', 'Stored residues')
plot(times, real(scaled_original_signals(selected_idx, :)), ':', ...
    'DisplayName', 'Stored residues, scaled')
plot(times, real(fitted_signals(selected_idx, :)), '-.', ...
    'DisplayName', 'Fitted residues')
xlabel('Time')
ylabel('Late-time signal')
grid on
legend('Location', 'best')

title(layout, sprintf('%s, receiver angle %.4f rad', ...
    object_label, selected_angle), 'Interpreter', 'none')
exportgraphics(fig, figure_path, 'ContentType', 'vector');
end


function make_angle_summary(angles, residue_error, scaled_residue_error, ...
        power_hellinger, original_signal_error, ...
        scaled_original_signal_error, fitted_signal_error, object_label, ...
        figure_path)
%MAKE_ANGLE_SUMMARY  Plot residue and reconstruction errors versus angle.

angle_degrees = mod(angles * 180 / pi, 360);
[angle_degrees, order] = sort(angle_degrees);

fig = figure('Color', 'w');
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile(layout)
semilogy(angle_degrees, max(residue_error(order), eps), ...
    'DisplayName', 'Raw residue error')
hold on
semilogy(angle_degrees, max(scaled_residue_error(order), eps), ...
    'DisplayName', 'Scale-aligned residue error')
semilogy(angle_degrees, max(power_hellinger(order), eps), ...
    'DisplayName', 'Relative-power Hellinger distance')
xlabel('Receiver angle (degrees)')
ylabel('Residue discrepancy')
xlim([0, 360])
grid on
legend('Location', 'best')

nexttile(layout)
semilogy(angle_degrees, max(original_signal_error(order), eps), ...
    'DisplayName', 'Stored residues')
hold on
semilogy(angle_degrees, max(scaled_original_signal_error(order), eps), ...
    'DisplayName', 'Stored residues, scaled')
semilogy(angle_degrees, max(fitted_signal_error(order), eps), ...
    'DisplayName', 'Fitted residues')
xlabel('Receiver angle (degrees)')
ylabel('Relative signal error')
xlim([0, 360])
grid on
legend('Location', 'best')

title(layout, object_label + " residue diagnostics", 'Interpreter', 'none')
exportgraphics(fig, figure_path, 'ContentType', 'vector');
end
