clearvars

%% Set up parameters
nsize = 2000;  % Number of neurons in the network
theta = 0;
tf = 4;
p = 200;       % Number of available stimulus pairs
b = 150;
mu = 0.97;

% Each row defines one condition for the first two stimulus pairs.
inputX = [1, 0; 0, 0; 0, 1; 1, 0];
inputY = [0, 0; 0, 1; 0, 1; 1, 0];

%% Run network simulations for match and mismatch conditions
imm = zeros(4, nsize);
x = zeros(p, 1);
y = zeros(p, 1);
h0 = zeros(nsize, 1);

mulnormal = mvnrnd([0, 0], [1, sqrt(mu); sqrt(mu), 1], nsize*p);
w = reshape(mulnormal(:, 1), nsize, p);
v = reshape(mulnormal(:, 2), nsize, p);
jmatrix = (w * w' + v * v') / nsize;

for i = 1:4
    x(:) = 0;
    y(:) = 0;
    x(1:2) = inputX(i, :)';
    y(1:2) = inputY(i, :)';

    dhdt = @(t, h) -h - b * jmatrix * max(h - theta, 0) + b * (w*x + v*y);
    [~, h] = ode45(dhdt, [0, tf], h0);
    imm(i, :) = h(end, :);
end

mmx1 = imm(1, :); % Mismatch response for 1st stimulus pair
m1 = imm(4, :);   % Match response for 1st stimulus pair
mmy2 = imm(2, :); % Mismatch response for 2nd stimulus pair
m2 = imm(3, :);   % Match response for 2nd stimulus pair

%% Save the data
outputFile = ['PairedRep_', '_b_', strtrim(num2str(b)), '_theta_', ...
    num2str(theta), '_p_', num2str(p), '.mat'];
save(outputFile, 'mu', 'theta', 'b', 'p', 'inputX', 'inputY', ...
    'mmx1', 'm1', 'mmy2', 'm2');
