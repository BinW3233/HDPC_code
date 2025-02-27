clear
%% Set up parameters
size = 2000;
theta = 0;
tf = 4; 
p = 200;
b = 150;
mu = 0.97;
inputX = [1,0; 0,0; 0,1; 1,0]; % Each row of the matrix gives the input x for thw two stimulus-pairs
inputY = [0,0; 0,1; 0,1; 1,0]; % Each row of the matrix gives the input y for thw two stimulus-pairs

%% Run network simulation for match and mismatch conditions for two stimulus-pairs (4 conditions in total)
imm = zeros(4, size);
x = zeros(p, 1);
y = zeros(p, 1);
h0 = zeros(1, size);
mulnormal = mvnrnd([0, 0], [1, sqrt(mu); sqrt(mu), 1], size*p);
w = reshape(mulnormal(:, 1), size, p);
v = reshape(mulnormal(:, 2), size, p);

for i = 1:4
    x(1:2) = inputX(i,:);
    y(1:2) = inputY(i,:);
    [t, h] = ode45(@(t, h) -h - b * (w * w' + v * v') * max(h - theta, 0)/size + b * (w * x + v* y), [0, tf], h0);
    imm(i, :) = h(end,:);
end

mmx1 = reshape(imm(1, :), [1, size]); % mismatch response for 1st stimulus-pair
m1 = reshape(imm(4, :), [1, size]);   % match response for 1st stimulus-pair
mmy2 = reshape(imm(2, :), [1, size]); % mismatch response for 2nd stimulus-pair
m2 = reshape(imm(3, :), [1, size]);   % match response for 2nd stimulus-pair

%% Save the data
title = ['PairedRep_', '_b_', strtrim(num2str(b)), '_theta_', num2str(theta), '_p_', num2str(p), '.mat'];
save(title,'mu','mmx1', 'mmy2', 'm2');