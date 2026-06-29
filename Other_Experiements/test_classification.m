noise_levels = [0];
% noise_levels = [0.0005];



% accuracies = [];
% for nl = noise_levels
%     nl
%     overall_acc = classify_objects("/scratch/msantana/time_domain_data_2D",0,nl,100,t_small, 250);
%     accuracies = [accuracies, overall_acc];
% end

%
% The center of the gaussian is 30 and it takes 11 seconds to turn off so
% we do the time after that.
%
%

accuracies = [];
for nl = noise_levels
    nl
    overall_acc = classify_objects("/home/msantana/Progamming/Target_Identification_2D/Data_Generation/plane_data",...,
        40,nl,359,[],[]);
    accuracies = [accuracies, overall_acc];
end

%%
accuracies = [];
for nl = noise_levels
    nl
    overall_acc = classify_objects("/home/msantana/Progamming/Target_Identification_2D/Developement_Tests/data_gen_planes_check",...,
        50,nl,359,[],[]);
    accuracies = [accuracies, overall_acc];
end

%% Is there actually a decrease in accuracy as I make the window length longer and longer?
noise_levels = [0.0007];
window_lengths = [20,50,80,100,150,200,250,300];
accuracies = [];
for wl = window_lengths
    for nl = noise_levels
        nl
        overall_acc = classify_objects("/scratch/msantana/time_domain_data_2D",0,nl,30,11.455033742029446, wl);
        accuracies = [accuracies, overall_acc];
    end
end

%%
figure(1)
clf
plot(window_lengths,accuracies,'-o')
xlabel("Window length (s)")
ylabel("Classification accuracy")
grid on



%%
ts = linspace(-20,20,1000);
w_0 = 50;
sigmas = 1.085736204758130;
inc_field_tol = 1e-16;

figure(1)
clf
plot(ts, real(at_gauss(ts,sigmas,w_0)));
hold on
t_small = gauss_picker_time(sigmas, inc_field_tol);
plot(t_small,real(at_gauss(t_small,sigmas,w_0)),'*');

%% How well does this actually work?
load("/scratch/msantana/time_domain_data_2D/circularcavity1.5.mat")
run_ind = 20;
figure(3)
clf
hold on
plot(ts,uff_all(20,:))
[~,peak] = max(abs(uff_all(20,:)));
[~,t_late_ind] = min(abs((ts(peak) + t_small) - ts));
t_late = ts(t_late_ind);
plot(t_late,uff_all(1,t_late_ind),'s','linewidth',5)
hold off

%%
function t_small = gauss_picker_time(sigmas, tol)
    t_small = sqrt(abs((4/sigmas) * log(2 * (sqrt(pi) / sqrt(sigmas)) * tol)));
end

function vals = at_gauss(t, sigmas, w_0)
    % This is the inverse fourier transform of exp(-(w - w_0)^2/sigmas)
    % sigmas is sigma squared
      vals = exp(-(1/4)*t .* (sigmas .* t + 4i * w_0));
      vals = vals / (sqrt(2) * sqrt(1/sigmas));
      vals = vals * (1/sqrt(2*pi));
end
