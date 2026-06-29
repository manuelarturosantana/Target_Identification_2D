%% Pole convergence seems less accurate than the scattered field convergence.
% It seems you should try to be really refined before capturing the
% poles...
%
%

figure(1)
clf
inds = [1,6];
for ii = length(inds)
    hold on
    plot(poles{inds(ii)},'*')
    hold off
end

%%
errs = [];
% poles = Poles;
conv_poles = poles{end};
for ii = 1:length(conv_poles)
    conv_pol = conv_poles(ii);
    loc_errs = [];
    for jj = 1:length(poles) - 1
        test_poles = poles{jj};
        loc_errs(end+1) = min(abs(test_poles - conv_pol));
    end
    errs = [errs; loc_errs];
end