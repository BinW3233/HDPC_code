clearvars

%% Parameters for the model
size = 600;
sigmaw = 1;
p = 2;
fw = 0.3;
theta = 0;
tf = 6;
mu = 0.97;
b = 150;
th = 7;
pL = @(x) (x .* (x >= 0 & x < th)) + (th .* (x >= th));

%% Specify the stimulus condition
x = zeros(p,1);
y = zeros(p,1);
x(1) = 1;
y(1) = 0;

%% Sample the random weight vectors
mulnormal = mvnrnd([0, 0], sigmaw^2*[1, sqrt(mu); sqrt(mu), 1], size*p);
w = reshape(mulnormal(:, 1), size, p);
v = reshape(mulnormal(:, 2), size, p);
jmatrix0 = (w * w' + v * v') / size;
rand_matrix = rand(size);
kij = rand_matrix < fw;
jmatrix = jmatrix0.*kij;

%% Get mismatch 1 response
h0 = zeros(size, 1);
dhdt = @(t, h) -h - b * jmatrix * pL(h - theta) + b * (w * x + v * y);
[t1, h1] = ode23tb(dhdt, [0, tf], h0);
r1 = pL(h1); % time-dependent firing rate

%% Save the data
title = 'Sparse_NN.mat';
save(title, 'h1', 'r1');


