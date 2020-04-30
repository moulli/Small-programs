%% 1-dimension

data = [randn(5,1)-10; 
        randn(5,1)+10];

phi = @(x) exp(-.5*x.^2)/sqrt(2*pi);

% true density
tpdf = @(x) phi(x+10)/2+phi(x-10)/2;

% kernel density estimation
bandwidth = std(data)*(4/3/numel(data))^(1/5);
kernel = @(x) mean(phi((x-data)/bandwidth)/bandwidth);
kpdf = @(x) arrayfun(kernel,x);

% plot
x = linspace(-25,+25,1000);
figure
hold on
plot(x, tpdf(x))
plot(x, kpdf(x))
legend('true density', 'kernel density estimation')


%% 2-dimension

mu1 = [-10, -10];
mu2 = [10, 10];
sigma = [1, 0; 0, 1];
data = [mvnrnd(mu1, sigma, 100);
        mvnrnd(mu2, sigma, 100)];
    
% true density
tpdf = @(x) mvnpdf(x, mu1, sigma)/2 + mvnpdf(x, mu2, sigma)/2;

% kernel density estimation
bandwidth = 3;
kernel = @(x) mean(mvnpdf((x-data)/bandwidth)/bandwidth);
kpdf = @(x) cellfun(kernel, mat2cell(x, ones(size(x, 1), 1)));

% plot
x = linspace(-25, +25, 200);
[X, Y] = meshgrid(x, x);
figure
hold on
plot3(X(:), Y(:), tpdf([X(:), Y(:)]))
plot3(X(:), Y(:), kpdf([X(:), Y(:)]))


%% 3-dimension

mu1 = [-10, -10, -10];
mu2 = [10, 10, 10];
sigma = [1, 0, 0; 0, 1, 0; 0, 0, 1];
data = [mvnrnd(mu1, sigma, 100);
        mvnrnd(mu2, sigma, 100)];
    
% true density
tpdf = @(x) mvnpdf(x, mu1, sigma)/2 + mvnpdf(x, mu2, sigma)/2;

% kernel density estimation
bandwidth = 3;
kernel = @(x) mean(mvnpdf((x-data)/bandwidth)/bandwidth);
kpdf = @(x) cellfun(kernel, mat2cell(x, ones(size(x, 1), 1)));

% plot
x = linspace(-25, +25, 75);
[X, Y, Z] = meshgrid(x, x, x);
figure
hold on
% T = tpdf([X(:), Y(:), Z(:)]);
% T = (max(T(:)) - T) ./ max(T(:));
% T_keep = (T < 0.999);
% scatter3(X(T_keep), Y(T_keep), Z(T_keep), [], T(T_keep).*(ones(1, 3)), '.')
K = kpdf([X(:), Y(:), Z(:)]);
K = (max(K(:)) - K) ./ max(K(:));
K_keep = (K < 0.75);
scatter3(X(K_keep), Y(K_keep), Z(K_keep), [], K(K_keep).*[0, 1, 1] + [1, 0, 0], '.')
axis equal




