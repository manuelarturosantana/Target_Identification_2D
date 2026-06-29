figure(1)
clf
rec_loc = 3;
inds = ts < 80;
plot(ts(inds),real(uff_all(rec_loc,inds)))

%%
figure(2)
plot(psp.ws,f_sols_all(rec_loc,:))

%%
figure(3)
ts2 = linspace(20,50,3000);
vals = psp.a_t(ts2);
figure(3)
clf
% semilogy(ts,abs(vals))
plot(ts2,real(vals))
title("Time Domain Gaussian Decay")
xlabel("ws")
ylabel("Magnitude of Gaussian")

%%
load('data_gen_planes_check/twinjet_plane4.mat')

inds = ts < 40;
for ii = 1:10
    figure(1)
    clf
    plot(ts(inds),real(uff_all(ii,inds)))
    pause(1)
end