clearvars

%% Define the model parameters
nsize = 600;
sigmaw = 1;
p = 2;
theta = 0;
mu = 0.97;
b = 150;
amp = 1;
tpre = 0.1;
txy = 10;
delta_t = 0.1; % inter-stimulus interval

%% Sample the random weight vectors
h0 = zeros(nsize, 1);
mulnormal = mvnrnd([0, 0], sigmaw^2*[1, sqrt(mu); sqrt(mu), 1], nsize*p);
w = reshape(mulnormal(:, 1), nsize, p);
v = reshape(mulnormal(:, 2), nsize, p);
jmatrix = (w * w' + v * v') / nsize;
x(1) = amp;
y(1) = amp;

%% 1st stage: no input
dhdt = @(t, h) (-h - b * jmatrix * max(h - theta, 0))/10; 
[t1, h1] = ode45(dhdt, [0, tpre], h0);

%% 2nd stage: pulse input x is presented initially

[t2, h2] = ode45(dhdt, [0, delta_t], h1(end,:)'+b*w*x);

%% 3rd stage: pulse input y is presented initially

[t3, h3] = ode45(dhdt, [0, txy], h2(end,:)'+b*v*y);

%% Combine time intervals
ht = cat(1,h1,h2,h3);
t_tot = cat(1, t1, t2+t1(end), t3+t2(end)+t1(end));

%% Voltage responses in the x only and y only conditions
[tx1, hx1] = ode45(dhdt, [0, tpre], h0);
[tx2, hx2] = ode45(dhdt, [0, delta_t + txy], hx1(end,:)'+b*w*x);

[ty1, hy1] = ode45(dhdt, [0, tpre+delta_t], h0);
[ty2, hy2] = ode45(dhdt, [0, txy], hy1(end,:)'+b*v*y);

tx = cat(1, tx1, tx2+tx1(end));
ty = cat(1, ty1, ty2+ty1(end));
hx = cat(1,hx1,hx2);
hy = cat(1,hy1,hy2);

%% save the data
title = ['TDcode_longtxy_','mu_', num2str(mu),'_b_', strtrim(num2str(b)), '_deltat_', num2str(delta_t), '_p_', num2str(p), '.mat'];
save(title, "t_tot", "ht", "tx", "hx", "ty", "hy", "jmatrix", "w", "v", "tpre", "delta_t", "txy", "amp", "b", "theta", "p", "mu");
