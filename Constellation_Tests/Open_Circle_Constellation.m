% Script to compute the eigenvalues of the open circle, and compare to the
% eigenvalues of various constellations

% Open Questions
% Is it really two poles, or is it one?
% What happens with two different scatterers
% Did we miss poles, or are they no longer poles

box_lims = [49.2, 49.6, -0.2, 0];

N = ceil((box_lims(2) * 12) / 2);

rp = rec_params("xlims",box_lims(1:2),"ylims",box_lims(3:4),"ntb",150,"nlr",150,"m_rec",7);
figure(1)
clf


%% Compute_Single Circle Eigenvalues
curve = Build_Curve("circular cavity", N, 0.75);

xshift = [];
yshift = [];
draw_constellation = true;
poles1 = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);


%%
xshift = [2.5];
yshift = [0];
poles2 = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);

%%
xshift = [2.5, 5.5];
yshift = [0, 0];
poles3 = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);
return
%%
% Next verify everything is a actually a true pole via the SVD. 
% Think about what this is going to mean for GLRT etc.
% save the data in an easy to present form for Oscar and Vicente

figure(1)
clf
hold on
plot(poles1,'*')
plot(poles2,'^')
plot(poles3,'p')

save("Poles123.mat","poles3","poles2","poles1")

%% Did we actually lose the poles, or did we miss some
[lp, mc] = setup_constellation(curve, xshift, yshift, draw_constellation)
%%
[~,I] = min(abs(49.4707 - poles1));
% [~,I] = min(abs(49.3735 - poles1));
test_p = poles1(I);


%% They are being missed
% I think this because the SVD shows two nearby poles (singular to 1e-8)
% and one right on top of it pole (singular to 1e-16)
[~,I] = mink(abs(test_p - poles3),1);
test_p_3 = poles3(I);
disp("SVD of Single Circle Pole")

S = svd(lp.bie_mat(test_p));
S(end)
format long
disp("SVDs of Three Circle Poles")
for ii = 1:length(test_p_3)
    S = svd(lp.bie_mat(test_p_3(ii)));
    S(end-3:end)
end

%%
test_zs = test_p + 0.01 * exp(2*pi*1i * (1:100)/100);
u = rand(1,N*3); v = randn(N*3,1);
f = @(k) u * (lp.bie_mat(k) \ v);

f_vals = [];
for k = test_zs
    f_vals(end+1) = f(k);
end

%%
[~, poles_near] = aaa(f_vals,test_zs,'tol',1e-7);

% figure(2)
% clf
% hold on
% plot(poles1(inpolygon(real(poles1),imag(poles1),real(test_zs),imag(test_zs))),'*')
% plot(poles2,'^')
% plot(poles3,'p')
% plot(poles_near,'s')

figure(2)
clf
hold on

idx1 = inpolygon(real(poles1), imag(poles1), real(test_zs), imag(test_zs));
idx2 = inpolygon(real(poles2), imag(poles2), real(test_zs), imag(test_zs));
idx3 = inpolygon(real(poles3), imag(poles3), real(test_zs), imag(test_zs));
idxn = inpolygon(real(poles_near), imag(poles_near), real(test_zs), imag(test_zs));

plot(poles1(idx1), '*')
plot(poles2(idx2), '^')
plot(poles3(idx3), 'p')
plot(poles_near(idxn), 's')

axis equal
grid on
hold off

%%
ref_ps = [];
for ii = 1:length(poles_near)
    p_ref = secant_method(@(k) 1 / f(k),poles_near(ii), poles_near(ii) + 1e-5, 10,1e-13,1);
    ref_ps(end+1) = p_ref;
    S = svd(lp.bie_mat(p_ref));
    S(end-3:end)
end

%%
p1 = ref_ps(3);
p2 = ref_ps(end-1);
abs(p1 - p2)

S = svd(lp.bie_mat(p1));
S(end)
S = svd(lp.bie_mat(p2));
S(end)
S = svd(lp.bie_mat(p1 + p2 / 2 ));
S(end)
%% It appears that this eigenvalue only gets split into two eigenvalues.






