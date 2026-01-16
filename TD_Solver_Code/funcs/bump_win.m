function vals = bump_win(t, H)
   % Window function as described by 12 in the paper by Anderson-Bruno-Lyon
   % Inputs:
   %      t : The values to be evaluated at
   %      H : The parameter H to be used in the evaluation.
   % Outputs:
   %    vals: The values of the field
   vals = t;
   % We need to check for cases where  there is not one of these cases.
   w = (t(t >= -H & t < -H/2) + H) / (H/2);
   if ~isempty(w)
        vals(t >= -H & t < -H/2) = 1 - eta(w); 
   end

   w = ones(size(t(t >= -H/2 & t <= H/2)));
   if ~isempty(w)
      vals(t >= -H/2 & t <= H/2) = w;           
   end

   w = (t(t >= H/2 & t < H) - H/2) / (H/2);
   if ~isempty(w) 
      vals(t >= H/2 & t < H)   = eta(w);     
   end

   w = 0 * vals(abs(t) >= H);
   
   if~isempty(w)
    vals(abs(t) >= H) = w;          
   end
end

function vals = eta(u)
    vals = 2 * exp(-1 ./ u) ./ (u - 1);
    vals = exp(vals);
end