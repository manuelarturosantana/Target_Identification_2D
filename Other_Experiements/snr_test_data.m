%%
% Some things which realy affect the size of the scattered field
% 1. The location. Being very much in the resonant region, vs being not
% makes a huge difference. In particular the exterior of the open circles
% is very much not trapping, and so nothing outside of the cone of the
% trapping region will see the resonant behavior very well.
% 2. Exciting enough poles.
% 3. Having higher frequency content.



%%

fw = 5;
gap = 1;
load("/scratch/msantana/snr_tests_data/circularcavity_fw" + num2str(fw) + "_gap" + num2str(gap) + ".mat")
for rec_ind = 1:10

    figure(rec_ind)
    clf
    tiledlayout(1,2)
    nexttile
    plot(curveX, curveY,'linewidth',2)
    hold on
    [theta,rho] = cart2pol(xs(rec_ind), ys(rec_ind));
    [x,y] = pol2cart(theta,rho*5);
    plot(x,y,'k*','markersize',15,'linewidth',5);
    legend("Scatterer","Reciever Location")
    axis equal
    
    nexttile
    inds = ts < 100 & ts > 20;
    sp_full = signal_power(ts(inds), uff_all(rec_ind,inds));
    inds_lt = ts > 40 & ts < 100;
    sp_lt   = signal_power(ts(inds_lt), uff_all(rec_ind, inds_lt));
    plot(ts(inds), real(uff_all(rec_ind,inds)))
    title("Full Signal Power " + num2str(sp_full) + " Late Time Signal Power " + num2str(sp_lt));
end



% [ind_start,ind_end] = signal_window_threshold(ts, uff_all(rec_ind,:), 5, 1e-6,1e-6);
% 
% 
% 
% 
% figure(2)
% clf
% plot(ts(ind_start:ind_end), real(uff_all(rec_ind,ind_start:ind_end)))




function P = signal_power(t, r)
    P = trapz(t, abs(r).^2);
end

function [ind_start, ind_end] = signal_window_threshold(t, r, window_size, tol_high, tol_low)
    % Find indices when first window exceeds tol_high, then when it drops below tol_low
    
    dt = t(2) - t(1);
    win_samples = round(window_size / dt);
    
    ind_start = [];
    ind_end = [];
    
    % Sliding window to find power
    for i = 1:length(t) - win_samples + 1
        window = r(i:i + win_samples - 1);
        P = trapz(t(i:i + win_samples - 1), abs(window).^2);
        
        if isempty(ind_start) && P > tol_high
            ind_start = i;
        end
        
        if ~isempty(ind_start) && isempty(ind_end) && P < tol_low
            ind_end = i;
            return
        end
    end
end
