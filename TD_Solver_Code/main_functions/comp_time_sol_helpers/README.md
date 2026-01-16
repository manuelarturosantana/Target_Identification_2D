## Comp_Tim_Sol Helpers
Unfortunately MATLAB's parfor loop makes using if statements inside of it really hairy. In
fact since everything has to be defined before the loop it is not doable. Therefore we have 
to split up the code and write many similar functions leading to code duplication (sadface).
This directory contains those codes.