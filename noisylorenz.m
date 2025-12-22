%% Lorenz attractor: time-series delay-coordinate embedding
%Bill Ostrie, Jac Ostrie, Ryan Caywood
clear; clc; close all;

function dXdt = lorenz_rhs(~, X)
    % Unpack state
    x = X(1);
    y = X(2);
    z = X(3);

    sigma = 10;
    beta = 8/3;
    rho = 28;

    % Lorenz equations
    X = sigma * (y - x);
    Y = x * (rho - z) - y;
    Z = x * y - beta * z;

    dXdt = [X; Y; Z];
end

%% 2. Time settings for integration
t0 = 0;          % initial time
tf = 60;         % final time
dt = 0.01;       % approximate time step for saving
tspan = t0:dt:tf;

%% 3. Initial condition and function handle (0,1,0), (1,1,1), (0,1.0,1.05)
x0 = 0;
y0 = 1;
z0 = 1.05;
X0 = [x0; y0; z0];

odefun = @(t,X) lorenz_rhs(t, X);
sol = ode45(odefun, [t0 tf], X0);

% Evaluate solution on tspan grid
X = deval(sol, tspan).';   % transpose to get [N x 3] matrix
t = tspan(:);              % column vector for times
x = X(:,1);
y = X(:,2);
z = X(:,3);

%% ADD NOISE
noise_level = 0.3;        % change for stronger/weaker noise
rng(1);                   % reproducible
X_noisy = x + noise_level * randn(size(x));
Y_noisy = y + noise_level * randn(size(y));
Z_noisy = z + noise_level * randn(size(z));

% original attractor in (x,y,z)
figure;
plot3(X_noisy, Y_noisy, Z_noisy, 'LineWidth', 0.5);
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('Original Noisy attractor in (x,y,z) space');
view(30, 20);
drawnow;





figure;
plot(X_noisy,'.r');
hold on;
%plot(y,'.g');
title('x coordinate time series');

figure;
plot(Y_noisy,'.r');
title('y coordinate time series');

writematrix(x,'Lorenz_attractor_x_noise');



% build delay coordinates
m = 3;
tau = 10;

N_points = size(X_noisy,1);
N_coords = N_points - m * tau;
delay_coords = zeros(m, N_coords);

for ii = 1:N_coords
    for jj = 1:m
        delay_coords(jj, ii) = X_noisy(ii + (m - jj) * tau);
    end
end

writematrix(delay_coords,'Noisy_Lorenz_delay_coords');

figure;
plot(delay_coords(1,:), '.b', 'MarkerSize', 3);
title('1D reconstruction');

% plot 2D reconstruction or projection
figure;
plot(delay_coords(1,:), delay_coords(2,:), '.r', 'MarkerSize', 3);
title('2D reconstruction');

% plot 3D reconstruction or projection
figure;
plot3(delay_coords(1,:), delay_coords(2,:), delay_coords(3,:), '.k', 'MarkerSize', 3);
title('3D reconstruction')



figure;
plot(delay_coords(1,:), '.b');
hold on;
plot(X_noisy,'.g');
title('first coordinate time series');



