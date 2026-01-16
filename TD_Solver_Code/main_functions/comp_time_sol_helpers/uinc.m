function val = uinc(x,y,kappa,t, w_0, sigmas)
    % Simple gaussian incident field following the convergence example in the paper.
    % sigmas is sigma squared
      r = dot(kappa, [x;y]);
      val = exp(-0.25 *(r - t) * (r * sigmas - sigmas * t - 4i * w_0));
      val = val / (sqrt(2) * sqrt(1/sigmas));
      % This came from wolfram Alpha, which by default uses the symmetric
      % normalization, so we follow the paper to have the inverse have the
      % 2pi normalization.
      val = val * (1/sqrt(2*pi));
end
