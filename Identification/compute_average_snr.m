function snr = compute_average_snr(data, model, params)
%COMPUTE_AVERAGE_SNR  Compute SNR statistics averaged across all objects.
%
%   snr = COMPUTE_AVERAGE_SNR(data, model, params) computes SNR for each
%   object at each receiver angle and noise level, then computes average
%   SNR and variance of SNR across all objects at each angle and noise level.
%
%   Output table contains:
%       sigma              noise standard deviation
%       receiver_angle     receiver angle
%       snr_full_db        SNR (full signal) per object
%       snr_late_db        SNR (late signal) per object
%       snr_full_mean_db   mean SNR (full signal) across objects
%       snr_late_mean_db   mean SNR (late signal) across objects
%       snr_full_var_db    variance of SNR (full signal) across objects
%       snr_late_var_db    variance of SNR (late signal) across objects

    nobj = numel(data);
    sigmas = params.noise_sigmas(:);
    nsigma = numel(sigmas);
    
    late_idx = model(1).late_idx;
    
    % Determine receiver angles from first object (assume all objects have same angles)
    if isempty(data(1).receiver_angles)
        receiver_angles = [1];
        nrec = 1;
    else
        receiver_angles = data(1).receiver_angles(:);
        nrec = numel(receiver_angles);
    end
    
    % Initialize 3D tensors: rows=sigma, columns=angle, depth=object
    snr_full_db_tensor = zeros(nsigma, nrec, nobj);
    snr_late_db_tensor = zeros(nsigma, nrec, nobj);
    snr_full_mean_db_tensor = zeros(nsigma, nrec);
    snr_late_mean_db_tensor = zeros(nsigma, nrec);
    snr_full_var_db_tensor = zeros(nsigma, nrec);
    snr_late_var_db_tensor = zeros(nsigma, nrec);
    
    % Loop over noise levels
    for isig = 1:nsigma
        sigma = sigmas(isig);
        
        % Loop over receiver angles
        for irec = 1:nrec
            
            snr_full_obj = zeros(nobj, 1);
            snr_late_obj = zeros(nobj, 1);
            
            % Loop over objects and compute SNR for each
            for jobj = 1:nobj
                f_full = data(jobj).signals(irec, :).';
                f_late = data(jobj).signals(irec, late_idx:end).';
                
                rms_full = sqrt(mean(abs(f_full).^2));
                rms_late = sqrt(mean(abs(f_late).^2));
                
                if sigma == 0
                    snr_full_obj(jobj) = Inf;
                    snr_late_obj(jobj) = Inf;
                else
                    snr_full_obj(jobj) = 20 * log10(rms_full / sigma);
                    snr_late_obj(jobj) = 20 * log10(rms_late / sigma);
                end
            end
            
            % Store SNR values in 3D tensor
            snr_full_db_tensor(isig, irec, :) = snr_full_obj;
            snr_late_db_tensor(isig, irec, :) = snr_late_obj;
            
            % Handle Inf values for mean/variance computation
            finite_full = isfinite(snr_full_obj);
            finite_late = isfinite(snr_late_obj);
            
            if any(finite_full)
                snr_full_mean_db_tensor(isig, irec) = mean(snr_full_obj(finite_full));
                snr_full_var_db_tensor(isig, irec) = var(snr_full_obj(finite_full));
            else
                snr_full_mean_db_tensor(isig, irec) = Inf;
                snr_full_var_db_tensor(isig, irec) = NaN;
            end
            
            if any(finite_late)
                snr_late_mean_db_tensor(isig, irec) = mean(snr_late_obj(finite_late));
                snr_late_var_db_tensor(isig, irec) = var(snr_late_obj(finite_late));
            else
                snr_late_mean_db_tensor(isig, irec) = Inf;
                snr_late_var_db_tensor(isig, irec) = NaN;
            end
        end
    end
    
    % Create output struct
    snr = struct( ...
        'sigma', sigmas, ...
        'receiver_angle', receiver_angles, ...
        'snr_full_db', snr_full_db_tensor, ...
        'snr_late_db', snr_late_db_tensor, ...
        'snr_full_mean_db', snr_full_mean_db_tensor, ...
        'snr_late_mean_db', snr_late_mean_db_tensor, ...
        'snr_full_var_db', snr_full_var_db_tensor, ...
        'snr_late_var_db', snr_late_var_db_tensor);
end
