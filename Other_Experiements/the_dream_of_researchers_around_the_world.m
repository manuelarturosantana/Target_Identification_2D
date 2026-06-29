load("/scratch/msantana/time_domain_data_2D/ellipticcavity1.5.mat")
fake_poles = pols;
%%
load("/scratch/msantana/time_domain_data_2D/circularcavity0.5.mat")

% figure(1)
uff_ind = 10;
% plot(ts,uff_all(uff_ind,:))

%%
nl = 1e-2;
signal = uff_all(uff_ind,:) + nl * randn(1, length(uff_all(uff_ind,:))) + 1i * nl * randn(1, length(uff_all(uff_ind,:))) ;
figure(1)
plot(ts,real(signal));
figure(2)
clf
hold on
plot(ts, log(abs(signal)))
[val,ind] = max(log(abs(signal)));
plot(ts(ind),abs(log(abs(signal(ind))),'*','markersize',10)
hold off
return






%%
estimated_start_time = 0;

decay_tol = 1e-17;
fastest_decay_rate = max(abs(imag(pols)));  % what if this pole is really weak...
window_size_sec = -log(decay_tol) / fastest_decay_rate;

signal = uff_all(uff_ind,:) + 10 * randn(1, length(uff_all(uff_ind,:))) + 1i * 10 * randn(1, length(uff_all(uff_ind,:))) ;
% figure(2)
% plot(ts,real(signal));

% start with the length of the signal that I took
ts_test = estimated_start_time:15:ts(end);

glrt_powers = [];
for ii = 1:length(ts_test)
    t_inds = (ts > ts_test(ii)) & (ts < (ts_test(ii) + window_size_sec));
    glrt_powers(end+1) = glrt(ts(t_inds), signal(t_inds), pols);
end

%%
fastest_decay_rate = max(abs(imag(fake_poles)));  % what if this pole is really weak...
window_size_sec = -log(decay_tol) / fastest_decay_rate;

% start with the length of the signal that I took
ts_test = estimated_start_time:15:ts(end);
glrt_powers_fake = [];
for ii = 1:length(ts_test)
    t_inds = (ts > ts_test(ii)) & (ts < (ts_test(ii) + window_size_sec));
    glrt_powers_fake(end+1) = glrt(ts(t_inds), signal(t_inds), fake_poles);
end

max(glrt_powers)
max(glrt_powers_fake)

%%
function result = glrt(tt, ft, poles)
% Generalized Likelihood Ratio method for object classification from poles

eps_thresh = 1e-12;  % drop small singular values

% Build Z matrix
Z = zeros(length(tt), length(poles));
for idx = 1:length(poles)
    Z(:, idx) = exp(-1j * poles(idx) * tt(:));
end

% Solve (Z^H Z)^(1/2) x = Z^H ft with the SVD
[U, S, ~] = svd(Z, 'econ');
s = diag(S);
keep = s > eps_thresh * max(s);
Ur = U(:, keep);

Uh_y = Ur' * ft(:);   % U* y  (.' is transpose, ' is conjugate transpose)
result = norm(Uh_y, 2) / length(tt);

end

%% AI Gen, not ensured that it is right.
% fs = 1000;
% t  = (0 : 1/fs : 1 - 1/fs);
% x  = chirp(t, 10, 1, 200);          % linear chirp 10 -> 200 Hz
% [W, t_ax, f_ax, tm, fm] = wigner_ville(x, fs);
% wigner_ville_plot(W, t_ax, f_ax, tm, fm);
function [W, t_axis, f_axis, t_marginal, f_marginal] = wigner_ville(signal, fs)
% WIGNER_VILLE  Compute the Wigner-Ville distribution of a time signal,
%               along with its time and frequency marginals.
%
% Inputs:
%   signal  - real or complex 1-D time-domain signal (column or row vector)
%   fs      - sampling frequency in Hz
%
% Outputs:
%   W           - Wigner-Ville distribution  [N x N], real-valued
%   t_axis      - time axis (seconds),       [1 x N]
%   f_axis      - frequency axis (Hz),       [1 x N]
%   t_marginal  - time marginal  = integral over f of W  [1 x N]
%   f_marginal  - frequency marginal = integral over t of W  [1 x N]
%
% Usage example:
%   fs = 1000;
%   t  = 0 : 1/fs : 1-1/fs;
%   x  = chirp(t, 0, 1, 200);          % linear chirp 0->200 Hz
%   [W, t_ax, f_ax, tm, fm] = wigner_ville(x, fs);
%   wigner_ville_plot(W, t_ax, f_ax, tm, fm);

    signal = signal(:);          % ensure column vector
    N      = length(signal);

    % --- analytic signal via Hilbert transform (suppresses cross-terms
    %     from negative frequencies for real inputs) ---
    if isreal(signal)
        z = hilbert(signal);
    else
        z = signal;
    end

    % --- time and frequency axes ---
    dt     = 1 / fs;
    t_axis = (0 : N-1) * dt;                      % [1 x N]
    f_axis = (0 : N-1) * (fs / N) - fs/2;         % centred, [1 x N]

    % --- WVD kernel: instantaneous autocorrelation ---
    % W(t, f) = integral_{-inf}^{inf} z(t + tau/2) * conj(z(t - tau/2))
    %           * exp(-j*2*pi*f*tau) dtau
    %
    % Discrete version: for each time index n, form the lag-product vector
    % R(n, m) = z(n+m) * conj(z(n-m))  over symmetric lags m = -(N-1)..N-1
    % then FFT over the lag index.

    W = zeros(N, N);    % rows = frequency, cols = time  (reshaped below)

    for n = 1 : N
        % symmetric lag range: keep only lags that stay inside the signal
        lags = 0 : min(n-1, N-n);   % non-negative lags (use symmetry)

        R = zeros(N, 1);             % zero-padded to length N for FFT

        for m = lags
            pos = n + m;             % z(t + tau/2) index  (1-based)
            neg = n - m;             % z(t - tau/2) index
            val = z(pos) * conj(z(neg));
            if m == 0
                R(1) = val;
            else
                % positive lag
                R(m + 1)     = val;
                % negative lag  (WVD kernel is Hermitian in lag)
                R(N - m + 1) = conj(val);
            end
        end

        % FFT over lag -> frequency slice; factor of 2 for one-sided kernel
        col        = 2 * real(fft(R));
        W(:, n)    = fftshift(col);  % centre zero-frequency
    end

    % W is [N_freq x N_time]; rows = frequency bins, cols = time samples
    % Both axes have length N.

    % --- marginals ---
    % time marginal   : integrate W over frequency  -> instantaneous power
    % frequency marginal: integrate W over time     -> power spectral density
    df          = fs / N;
    t_marginal  = sum(W, 1) * df;          % [1 x N]
    f_marginal  = sum(W, 2)' * dt;         % [1 x N]

end


% -------------------------------------------------------------------------
function wigner_ville_plot(W, t_axis, f_axis, t_marginal, f_marginal)
% Convenience plotting function for the WVD and its marginals.

    figure('Name', 'Wigner-Ville Distribution', 'NumberTitle', 'off', ...
           'Position', [100 100 900 700]);

    % --- main WVD panel ---
    ax_main = subplot(3, 3, [1 2 4 5]);
    imagesc(t_axis, f_axis, W);
    axis xy;
    colormap(ax_main, 'jet');
    colorbar;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Wigner-Ville Distribution');

    % --- time marginal (bottom) ---
    ax_t = subplot(3, 3, [7 8]);
    plot(t_axis, t_marginal, 'b', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Power');
    title('Time Marginal  \int W\,df');
    xlim([t_axis(1) t_axis(end)]);
    grid on;

    % --- frequency marginal (right) ---
    ax_f = subplot(3, 3, [3 6]);
    plot(f_marginal, f_axis, 'r', 'LineWidth', 1.2);
    ylabel('Frequency (Hz)');
    xlabel('Power');
    title('Freq. Marginal  \int W\,dt');
    ylim([f_axis(1) f_axis(end)]);
    grid on;

    % link axes for zoom/pan consistency
    linkaxes([ax_main ax_t], 'x');
    linkaxes([ax_main ax_f], 'y');

end