clearvars

%% Parameters for the model
nsize = 600;   % Number of neurons in the network
sigmaw = 1;
p = 2;         % Number of stimulus pairs
fw = 0.3;      % Connection probability for the sparse recurrent matrix
theta = 0;
tf = 6;
mu = 0.97;
b = 150;
th = 7;
pL = @(x) (x .* (x >= 0 & x < th)) + (th .* (x >= th));

%% Specify the stimulus condition
x = zeros(p, 1);
y = zeros(p, 1);
x(1) = 1;
y(1) = 0;

%% Sample the random weight vectors
mulnormal = mvnrnd([0, 0], sigmaw^2*[1, sqrt(mu); sqrt(mu), 1], nsize*p);
w = reshape(mulnormal(:, 1), nsize, p);
v = reshape(mulnormal(:, 2), nsize, p);
jmatrix0 = (w * w' + v * v') / nsize;
rand_matrix = rand(nsize);
kij = rand_matrix < fw;
jmatrix = jmatrix0 .* kij;

%% Get mismatch 1 response
h0 = zeros(nsize, 1);
dhdt = @(t, h) -h - b * jmatrix * pL(h - theta) + b * (w*x + v*y);
[t1, h1] = ode23tb(dhdt, [0, tf], h0);
r1 = pL(h1 - theta); % Time-dependent firing rate

%% Save the data
outputFile = 'Sparse_NN.mat';
save(outputFile, 't1', 'h1', 'r1', 'jmatrix', 'w', 'v', ...
    'theta', 'tf', 'mu', 'b', 'p', 'fw');
