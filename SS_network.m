clearvars
tic;

%% Set up parameters 
size = 500; % Number of neurons in the network
sigmaw = 1; 
p = 2;      % Number of stimulus-pair
theta = 0;  % Threshold for neural firing
tf = 5;     % Total simulation time
mu = 0.7;   % Learning stage
b = 150;    % Parameter b in the model  
h0 = zeros(size, 1);
x = zeros(p,1);
y = zeros(p,1);

%% Stimulus input
x(1) = 1;   % Strength of stimulus-pair input
y(1) = 0;

%% Generate random weight vectors
mulnormal = mvnrnd([0, 0], sigmaw^2*[1, sqrt(mu); sqrt(mu), 1], size*p);
w = reshape(mulnormal(:, 1), size, p);
v = reshape(mulnormal(:, 2), size, p);
jmatrix = (w * w' + v * v') / size;     % Recurrent weight matrix

%% Define the differential equation for neural dynamics and solve using ode45
dhdt = @(t, h) -h - b * jmatrix * max(h - theta, 0) + b * (w * x + v * y);
[t, h] = ode45(dhdt, [0, tf], h0);

%% Compute steady state voltage level and firing rate
hf = h(end,:);              % Steady state voltage levels
rf = max(hf - theta, 0); 
%plot(t, mean(h,2)) %plot the data

%% Save the data
title = 'ss_response.mat';
save(title, 'hf', 'rf');