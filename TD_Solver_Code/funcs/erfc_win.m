function vals = erfc_win(t, H)
   % Window function based on the complementary error function.
   % Inputs:
   %      t : The values to be evaluated at
   %      H : The parameter H to be used in the evaluation corresponding to window size.
   % Outputs:
   %    vals: The values of the field

   rho   = 5.805; % Leads to erfc_wind(t,t) = 1.1x10^-16.
   alpha = 0.5;   % The function will be one in |z| < \alpha H

   vals = t;

   vals(abs(t) < alpha * H) = 1;
   vals(abs(t) > H) = 0;

   erfc_inds = abs(t) >= alpha * H & abs(t) <= H;
   erfc_arg  = (abs(t(erfc_inds)) - alpha * H) / ((1 - alpha) * H);

   vals(erfc_inds) = (1/2) * erfc(-rho + 2 * rho * erfc_arg);

end