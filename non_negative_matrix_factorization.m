% Set focus
addpath(genpath('/home/ljp/Programs'))
root = '/home/ljp/Science/GeoffreysComputer/Projects/RLS/';
study = '';
date = '2019-10-01';
run = 06;
F = NT.Focus(root, study, date, run);

% Load Mmap
temp = adapted4DMatrix(F, 'corrected');

% Take layer, cast and clean Mmap
t = temp(:, :, 12, :);
t = permute(t, [1, 2, 4, 3]);
t = double(t);
t = t(25:end, 25:end, :);

% Take background and compute mean trace
background = t(300:end, 1:180, :);
t_back = mean(background);
t_back = mean(permute(t_back, [2, 3, 1]));
t_back = permute(t_back, [2, 1]);

% Linear regression on background to get temporal background matrix
tempo = fitlm((1:3000)', t_back');
ab = tempo.Coefficients.Estimate(1).*ones(size(t_back'))+tempo.Coefficients.Estimate(2).*(1:3000);

% Get subpart of t
npix = 30;
startx = 121;
starty = 271;
t = t(startx:startx+npix-1, starty:starty+npix-1, :);
% Plot subpart
figure
image(t(:, :, 1), 'CDataMapping', 'scaled')
axis equal

% Rearranging t to get F
F = zeros(size(t, 1)*size(t, 2), size(t, 3));
for i = 1:size(t, 3); ttemp = t(:, :, i); F(:, i) = ttemp(:); end

% Initialize S and A
numneu = npix;
S = rand(size(F, 1), numneu) * 10;
A = rand(numneu, size(F, 2)) * 10;

% Launch algorithm
rmse = sqrt(sum(sum((F - S*A).^2)));
fprintf('Initial rmse is %.3f \n', rmse);
for i = 1:100
    sb = (ab * ab') \ (F - S * A) * ab';
    sb(sb < 0) = 0.01;
    A = (S' * S) \ S' * (F - sb * ab);
    A(A < 0) = 0.01;
    S = ((A * A') \ A * (F - sb * ab)')';
    S(S < 0) = 0.01;
    S = (S - mean(S)) ./ std(S);
    rmse = sqrt(sum(sum((F - S*A - sb*ab).^2)));
    fprintf('Iteration %.0f, rmse is %.3f \n', [i, rmse]);
end

% Plot one so called neuron
Sn = zeros(npix, npix, numneu);
for i = 1:numneu
    Sn(:, :, i) = reshape(S(:, i), npix, npix);
end
figure
image(Sn(:, :, 1), 'CDataMapping', 'scaled')
axis equal
Sni = Sn(:, :, 1);
figure
hist(Sni(:), 50)
figure
image(mean(Sn, 3), 'CDataMapping', 'scaled')
axis equal
Sni = Sn;
Sni(Sni < 0) = 0;
figure
image(mean(Sni, 3), 'CDataMapping', 'scaled')
axis equal


