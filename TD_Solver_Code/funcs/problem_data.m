function ps = problem_data(varargin)
    % wrapper function which generates the problem data. Same arguments as
    % gen params.
    ps = gen_params(varargin{:});

    ps = init_params(ps);
end