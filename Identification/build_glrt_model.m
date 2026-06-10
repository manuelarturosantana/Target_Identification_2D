function model = build_glrt_model(data, params)
%BUILD_ID_MODEL  Build GLRT subspace model from pole data.
%
%   model = BUILD_ID_MODEL(data, params) builds one GLRT basis per target
%   configuration. The basis is obtained from the columns of
%
%       Z(t, p) = exp(-1i * t * p),
%
%   restricted to the late-time window.

    nobj = numel(data);

    if nobj == 0
        error('build_id_model:EmptyData', ...
              'Input data struct is empty.');
    end

    ts = data(1).ts(:);

    late_idx = find(ts >= params.late_time, 1, 'first');

    if isempty(late_idx)
        error('build_id_model:InvalidLateTime', ...
              'params.late_time is larger than the final time sample.');
    end

    tt = ts(late_idx:end);

    model = struct([]);

    for j = 1:nobj
        poles = data(j).poles(:).';

        Z = exp(-1i * tt(:) * poles);

        [U, S, ~] = svd(Z, 'econ');
        s = diag(S);

        if isempty(s)
            keep = false(size(s));
        else
            keep = s > params.svd_tol * max(s);
        end

        model(j).name = data(j).name;
        model(j).file_name = data(j).file_name;
        model(j).file_path = data(j).file_path;
        model(j).poles = data(j).poles(:);
        model(j).late_idx = late_idx;
        model(j).time_window = tt(:);
        model(j).Ur = U(:, keep);
        model(j).singular_values = s;
        model(j).retained_rank = nnz(keep);
    end

end
