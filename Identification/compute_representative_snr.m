function snr = compute_representative_snr(data, model, params)

    % Choose representative object.
    if strlength(string(params.snr_object_name)) == 0
        obj_idx = 1;
    else
        names = strings(numel(data), 1);

        for j = 1:numel(data)
            names(j) = string(data(j).name);
        end

        obj_idx = find(names == string(params.snr_object_name), 1);

        if isempty(obj_idx)
            error('classify_glrt_data:InvalidSnrObjectName', ...
                  'Could not find SNR object name "%s".', ...
                  string(params.snr_object_name));
        end
    end

    % Choose representative receiver.
    if isempty(params.snr_receiver_angle)
        rec_idx = 1;
    else
        if isempty(data(obj_idx).receiver_angles)
            error('compute_representative_snr:MissingReceiverAngles', ...
                  ['Cannot select SNR receiver by angle because ' ...
                   'receiver_angles is missing for object "%s".'], ...
                  data(obj_idx).name);
        end
    
        target_angle = params.snr_receiver_angle;
        receiver_angles = data(obj_idx).receiver_angles(:);
    
        angle_dist = abs(angle(exp(1i * (receiver_angles - target_angle))));
        [~, rec_idx] = min(angle_dist);
    end

    late_idx = model(1).late_idx;
    sigmas = params.noise_sigmas(:);
    nsigma = numel(sigmas);

    f_full = data(obj_idx).signals(rec_idx, :).';
    f_late = data(obj_idx).signals(rec_idx, late_idx:end).';

    rms_full = sqrt(mean(abs(f_full).^2));
    rms_late = sqrt(mean(abs(f_late).^2));

    snr_full_db = zeros(nsigma, 1);
    snr_late_db = zeros(nsigma, 1);

    for k = 1:nsigma
        sigma = sigmas(k);

        if sigma == 0
            snr_full_db(k) = Inf;
            snr_late_db(k) = Inf;
        else
            snr_full_db(k) = 20 * log10(rms_full / sigma);
            snr_late_db(k) = 20 * log10(rms_late / sigma);
        end
    end
    
    if isempty(data(obj_idx).receiver_angles)
        selected_receiver_angle = NaN;
    else
        selected_receiver_angle = data(obj_idx).receiver_angles(rec_idx);
    end

    snr = table( ...
        sigmas, ...
        snr_full_db, ...
        snr_late_db, ...
        repmat(obj_idx, nsigma, 1), ...
        repmat(rec_idx, nsigma, 1), ...
        repmat(string(data(obj_idx).name), nsigma, 1), ...
        repmat(selected_receiver_angle, nsigma, 1), ...
        'VariableNames', { ...
            'sigma', ...
            'snr_full_db', ...
            'snr_late_db', ...
            'object_idx', ...
            'receiver_idx', ...
            'object_name', ...
            'receiver_angle'});
end