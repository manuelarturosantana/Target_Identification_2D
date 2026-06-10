function overall_acc = classify_objects(data_dir, t_start, noise_level, n_signals, peak_offset)
% CLASSIFY_OBJECTS  Load .mat files, classify each signal row using GLRT,
%                   and report overall accuracy.
%
%   classify_objects(data_dir, t_start, noise_level)
%   classify_objects(data_dir, t_start, noise_level, n_signals)
%   classify_objects(data_dir, t_start, noise_level, n_signals, peak_offset)
%
%   Required inputs:
%     data_dir     - path to folder containing .mat files
%                    each file must have: pols, ts, uff_all
%     t_start      - only use time samples t >= t_start for classification
%                    (used as a fixed floor, or ignored per-signal when
%                    peak_offset is supplied)
%     noise_level  - standard deviation of additive complex Gaussian noise
%                    (e.g. 1e-3)
%
%   Optional inputs:
%     n_signals    - classify only the first n_signals rows of uff_all in
%                    each file. Omit or pass [] to use all rows.
%     peak_offset  - when supplied, the start time for each signal is
%                    computed as:
%                        t_start_i = t_peak_i + peak_offset
%                    where t_peak_i is the time of the absolute-value peak
%                    of that signal. t_start is still used as a lower bound:
%                        t_start_i = max(t_start_i, t_start)
%                    Omit or pass [] to use the fixed t_start for all signals.
%
%   Each .mat file represents one "true" object class.
%   Classification assigns a signal to the class whose poles yield the
%   highest GLRT statistic.

% -------------------------------------------------------------------------
% 0.  Parse optional arguments
% -------------------------------------------------------------------------
if nargin < 4 || isempty(n_signals)
    n_signals = Inf;   % use all rows
end
if nargin < 5 || isempty(peak_offset)
    peak_offset = [];  % fixed t_start mode
end
use_peak_mode = ~isempty(peak_offset);

% -------------------------------------------------------------------------
% 1.  Discover all .mat files
% -------------------------------------------------------------------------
fprintf('Settings: noise_level=%.2g | t_start=%.4g | n_signals=%s | peak_mode=%s\n\n', ...
        noise_level, t_start, ...
        num2str(n_signals), mat2str(use_peak_mode));
listing = dir(fullfile(data_dir, '*.mat'));
if isempty(listing)
    error('No .mat files found in: %s', data_dir);
end
num_classes = numel(listing);
fprintf('Found %d .mat file(s) in "%s".\n\n', num_classes, data_dir);

% -------------------------------------------------------------------------
% 2.  First pass – load and store the poles for every class
% -------------------------------------------------------------------------
all_poles = cell(num_classes, 1);
class_names = cell(num_classes, 1);

for k = 1:num_classes
    fpath = fullfile(listing(k).folder, listing(k).name);
    S = load(fpath, 'pols');
    all_poles{k} = S.pols(:);          % store as column vector
    class_names{k} = listing(k).name;
end

% -------------------------------------------------------------------------
% 3.  Second pass – classify every signal row in every file
% -------------------------------------------------------------------------
total_signals  = 0;
total_correct  = 0;

for true_class = 1:num_classes
    fpath = fullfile(listing(true_class).folder, listing(true_class).name);
    S = load(fpath, 'ts', 'uff_all');
    ts      = S.ts(:)';          % 1 x num_t  (row)
    uff_all = S.uff_all;         % n x num_t

    % --- cap to first n_signals rows ------------------------------------
    n = min(size(uff_all, 1), n_signals);
    class_correct = 0;

    for row = 1:n
        clean_full = uff_all(row, :).';   % full-length column vector

        % --- determine this signal's start time -------------------------
        if use_peak_mode
            [~, peak_idx] = max(abs(clean_full));
            t_peak        = ts(peak_idx);
            t_start_row   = max(t_peak + peak_offset, t_start);
        else
            t_start_row = t_start;
        end

        % --- select the time window t >= t_start_row --------------------
        time_mask = ts >= t_start_row;
        if ~any(time_mask)
            warning('File "%s", row %d: no samples >= t_start=%.4g. Skipping.', ...
                    class_names{true_class}, row, t_start_row);
            continue
        end
        tt           = ts(time_mask)';         % column vector of selected times
        clean_signal = clean_full(time_mask);  % already a column
        noise = (noise_level / sqrt(2)) * ...
                (randn(size(clean_signal)) + 1j * randn(size(clean_signal)));
        noisy_signal = clean_signal + noise;

        % --- evaluate GLRT against every class's poles ------------------
        scores = zeros(num_classes, 1);
        for c = 1:num_classes
            scores(c) = glrt(tt, noisy_signal, all_poles{c});
        end

        % --- classify as the class with the highest score ---------------
        [~, predicted_class] = max(scores);

        if predicted_class == true_class
            class_correct = class_correct + 1;
        end
    end

    % --- per-file report ------------------------------------------------
    fprintf('True class %2d (%s):\n', true_class, class_names{true_class});
    fprintf('  Signals: %d  |  Correct: %d  |  Accuracy: %.1f%%\n\n', ...
            n, class_correct, 100 * class_correct / n);

    total_signals = total_signals + n;
    total_correct = total_correct + class_correct;
end

% -------------------------------------------------------------------------
% 4.  Overall accuracy
% -------------------------------------------------------------------------
overall_acc = 0;
if total_signals > 0
    overall_acc = 100 * total_correct / total_signals;
    fprintf('==============================================\n');
    fprintf('Overall accuracy: %d / %d  (%.2f%%)\n', ...
            total_correct, total_signals, overall_acc);
    fprintf('==============================================\n');
else
    fprintf('No signals were classified (check t_start value).\n');
end

end % function classify_objects


% =========================================================================
%  Helper: GLRT statistic
% =========================================================================
function result = glrt(tt, ft, poles)
% GLRT  Generalised Likelihood Ratio statistic for one signal and one
%       set of poles.
%
%   result = glrt(tt, ft, poles)
%
%   tt    - column vector of time samples
%   ft    - column vector of signal values (same length as tt)
%   poles - column vector of complex poles

eps_thresh = 1e-12;   % drop singular values smaller than this fraction of max

% Build the Vandermonde-like Z matrix: Z(i,j) = exp(-1j * poles(j) * tt(i))
Z = zeros(length(tt), length(poles));
for idx = 1:length(poles)
    Z(:, idx) = exp(-1j * poles(idx) * tt(:));
end

% Economy SVD of Z
[U, S, ~] = svd(Z, 'econ');
s    = diag(S);
keep = s > eps_thresh * max(s);
Ur   = U(:, keep);           % retain only significant left singular vectors

% Project signal onto the column space of Z
Uh_y   = Ur' * ft(:);        % conjugate transpose projection
result = norm(Uh_y, 2);      % GLRT score
end