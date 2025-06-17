clearvars

%% Parameters for the model
size = 400;
sigmaw = 1;
p = 2;
theta = 0;
mu = 0.96;
b1 = 50; % parameter b for M1 and M3
b2 = 190; % parameter b for M2

%% Define the stimulus condition
x = zeros(p,1);
y = zeros(p,1);
x(1) = 1;
y(1) = 1;

%% Sample the random weight vectors
h0 = zeros(3*size, 1);
Cor0 = (1-sqrt(mu))*eye(10);
Cor1 = sqrt(mu)*ones(10);
mulnormal = mvnrnd(zeros(1,10), sigmaw^2*(Cor0+Cor1), size*p);

w1 = reshape(mulnormal(:, 1), size, p);
v3 = reshape(mulnormal(:, 2), size, p);
w2l = reshape(mulnormal(:, 3), size, p);
w2r = reshape(mulnormal(:, 4), size, p);
w3l = reshape(mulnormal(:, 5), size, p);
w3r = reshape(mulnormal(:, 6), size, p);
v1l = reshape(mulnormal(:, 7), size, p);
v1r = reshape(mulnormal(:, 8), size, p);
v2l = reshape(mulnormal(:, 9), size, p);
v2r = reshape(mulnormal(:, 10), size, p);


w2 = w2l * w2l'/size;
w3 = w3l * w3l'/size;
v1 = v1l * v1l'/size;
v2 = v2l * v2l'/size;

jmatrix1 = (w1 * w1'  + v1l *  v1l')/size;
jmatrix2 = (w2l *  w2l' + v2l *  v2l')/size;
jmatrix3 = (w3l *  w3l'+ v3 * v3')/size;

%% Define the differential equation and solve using ode45
dhdt = @(t, h) odefun(t, h, size, jmatrix1, jmatrix2, jmatrix3, w1, w2, w3, v1, v2, v3, b1, b2, theta, x, y);

[t, h] = ode23tb(dhdt, [0, 4], h0);
hf = h(end,:)';
rf = max(hf - theta, 0);

%% Save the data
title = 'Hierachical_response.mat';
save(title, 'hf', 'rf');


%% Define the neural dynamics for hierachical network
function dydt = odefun(t, h, size, jmatrix1, jmatrix2, jmatrix3, w1, w2, w3, v1, v2, v3, b1, b2, theta, x,y)
dydt = zeros(3*size,1);
dydt(1:size) = -h(1:size) - b1*jmatrix1*max(h(1:size)- theta, 0) + b1*(w1*x+v1*max(h(size+1:2*size)-theta,0));
dydt(size+1:2*size) = -h(size+1:2*size) - b2*jmatrix2*max(h(size+1:2*size)- theta, 0) + b2*(w2*max(h(1:size)-theta,0) + v2*max(h(2*size+1:3*size)-theta,0));
dydt(2*size+1:3*size) = -h(2*size+1:3*size) - b1*jmatrix3*max(h(2*size+1:3*size)- theta, 0) + b1*(w3*max(h(size+1:2*size)-theta,0)+v3*y);
end
