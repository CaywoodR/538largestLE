%% Lorenz attractor: time-series with additive Gaussian noise
clear; clc; close all;

function dXdt = lorenz_rhs(~, X)
    sigma = 10;
    beta  = 8/3;
    rho   = 28;

    x = X(1); y = X(2); z = X(3);
    dXdt = [sigma * (y - x);
            x * (rho - z) - y;
            x * y - beta * z];
end

%% ODE45
t0 = 0; tf = 1000; dt = 0.01;
tspan = t0:dt:tf;
X0 = [0; 1; 1.05];

odefun = @(t,X) lorenz_rhs(t, X);
sol = ode45(odefun, [t0 tf], X0);
X_full = deval(sol, tspan).';   % N×3 matrix
t = tspan(:);

%% Extract
X = X_full(:,1);
Y = X_full(:,2);
Z = X_full(:,3);

%% ADD NOISE
noise_level = 0.3;        % change for stronger/weaker noise
rng(1);                   % reproducible
X_noisy = X + noise_level * randn(size(X));
Y_noisy = Y + noise_level * randn(size(Y));
Z_noisy = Z + noise_level * randn(size(Z));

%% Plots
figure;
subplot(1,2,1)
plot3(X, Y, Z, 'k', 'LineWidth', 0.4)
title('Original Lorenz Attractor')
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on; view(30,20)

subplot(1,2,2)
plot3(X_noisy, Y_noisy, Z_noisy, 'r', 'LineWidth', 0.4)
title('Noisy Lorenz Attractor')
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on; view(30,20)

%% Output Noisy
LorenzData.t = t;
LorenzData.X = X_noisy;
LorenzData.Y = Y_noisy;
LorenzData.Z = Z_noisy;

%% ---------- Delay-Coordinate Embedding (Takens) ----------
% Goal: reconstruct state space from a single noisy observable (X_noisy)
% using an m-dimensional embedding with delay tau.

% ---- User choices ----
m = 3.0;                 % embedding dimension
tau_seconds = 0.55;    % delay in SECONDS (set this)
% If you want tau in samples instead, comment the line above and set:
 %tau_samples = 10;

% ---- Convert delay to samples (if given in seconds) ----
if exist('tau_samples','var') ~= 1
    % dt was defined earlier (sampling interval for tspan)
    tau_samples = max(1, round(tau_seconds / dt));
end

% ---- Build embedding from the scalar time series X_noisy ----
x = X_noisy(:);
x2 =X(:);% ensure column vector
N = numel(x2);                            % total samples
N_coords = N - (m-1)*tau_samples;        % number of embedded points

if N_coords < 1
    error('Embedding parameters too large: decrease tau or m.');
end

% Delay matrix: rows = points, columns = coordinates
delay_coords = zeros(N_coords, m);
for k = 1:m
    idx_start = 1 + (k-1)*tau_samples;
    idx_end   = idx_start + N_coords - 1;
    delay_coords(:, k) = x2(idx_start:idx_end);
end
t_emb = t(1:N_coords); 
% ---- Plots ----
figure;
plot(t,X_noisy)

figure

% -------------------------------------------------
% Subplot 1
% -------------------------------------------------
subplot(3,1,1)
plot(t, X, 'LineWidth', 1.5)
xlabel('t')
xlim([0 50])
ylabel('X')
title('X vs t')
grid on

% -------------------------------------------------
% Subplot 2
% -------------------------------------------------
subplot(3,1,2)
plot(t, Y, 'LineWidth', 1.5)
xlabel('t + 3')
ylabel('X')
xlim([0 50])
title('X vs (t + 3)')
grid on

% -------------------------------------------------
% Subplot 3
% -------------------------------------------------
subplot(3,1,3)
plot(t_emb, delay_coords(:,3), 'LineWidth', 1.5)
xlabel('t + 5')
xlim([0 50])
ylabel('x')
title('x vs (t + 5)')
grid on

figure;
plot(delay_coords(:,1), delay_coords(:,2), '.', 'MarkerSize', 3);
grid on; axis tight
xlabel('X(t)')
ylabel(sprintf('X(t + %d\\,samples)', tau_samples))
title(sprintf('2D Delay Reconstruction (m=%d, \\tau=%d samples)', m, tau_samples))


disp('✅ Delay embedding complete: "delay_coords" is [N_coords × m].');

%% ---------- Better autocorrelation + automatic tau candidates ----------
x = X_noisy(:);
x = x - mean(x);

max_lag_seconds = 5;
max_lag_samples = round(max_lag_seconds / dt);

[acf_vals, lags] = autocorr(x, 'NumLags', max_lag_samples);
lags_seconds = lags * dt;

% --- Candidate 1: first lag where |acf| <= 1/e (ignore lag 0) ---
thr = 1/exp(1);
idx_1e = find(abs(acf_vals(2:end)) <= thr, 1, 'first') + 1;  % +1 for lag indexing
tau_1e_samples = lags(idx_1e);
tau_1e_seconds = tau_1e_samples * dt;

% --- Candidate 2: first local minimum after lag 0 ---
mins = islocalmin(acf_vals);
mins(1) = false; % ignore lag 0
idx_min = find(mins, 1, 'first');
if isempty(idx_min)
    tau_min_samples = NaN;
    tau_min_seconds = NaN;
else
    tau_min_samples = lags(idx_min);
    tau_min_seconds = tau_min_samples * dt;
end

% --- Plot ACF with markers ---
figure;
plot(lags_seconds, acf_vals, 'LineWidth', 1.2); hold on;
yline(0,'--');
yline(thr,'--', '1/e');
yline(-thr,'--', '-1/e');
grid on; axis tight;
xlabel('\tau (seconds)');
ylabel('ACF');
title(sprintf('Autocorrelation of X (0 to %.1f s)', max_lag_seconds));



%% -------------------------------------------------------
%  Simple visual validation of tau (2D delay plots)
% --------------------------------------------------------

x = X_noisy(:);        % ensure column
x = x - mean(x);       % demean (important but simple)

taus = [1 2 4 5 10 40];   % candidate delays (in samples)

figure;

for i = 1:length(taus)
    tau = taus(i);

    % build delayed pairs
    x1 = x(1+tau:end);
    x2 = x(1:end-tau);

    subplot(2,3,i)
    plot(x1, x2, '.', 'MarkerSize', 3)
    axis equal tight
    grid on

    xlabel('x_n')
    ylabel(sprintf('x_{n-%d}', tau))
    title(sprintf('\\tau = %d', tau))
end

sgtitle('Visual validation of delay \tau via 2D delay plots')



%% ---------- Box-counting dimension (3D Lorenz or delay embedding) ----------

D_lorenz = [X_noisy, Y_noisy, Z_noisy];

% ---- Pick one ----
D = D_lorenz;   % <-- change to P_embed if desired

% ---- (1) Remove transient & (optional) downsample to speed up ----
t_transient = 10;                          % seconds to discard
idx0 = find(t >= t_transient, 1, 'first'); % first index after transient
D = D(idx0:end, :);

ds = 1;                 % downsample factor (increase if slow)
D = D(1:ds:end, :);

% ---- (2) Compute box counts over epsilons ----
eps_list = logspace(-4, 0.2, 40);  % adjust range if needed
% Note: eps range depends on scaling of D; we'll normalize below.

% Normalize to a unit-ish bounding box scale to make eps meaningful
Dmin = min(D, [], 1);
Dmax = max(D, [], 1);
Drng = Dmax - Dmin;
Dnorm = (D - Dmin) ./ Drng;        % now roughly inside [0,1]^d

N_eps = zeros(size(eps_list));



for i = 1:numel(eps_list)
    eps = eps_list(i);

    % Grid coordinates (box indices)
    idx = floor(Dnorm / eps);

    % Count occupied boxes (unique rows)
    N_eps(i) = size(unique(idx, 'rows'), 1);
end

% ---- (3) Fit slope over a chosen scaling window ----
xfit = log(1 ./ eps_list(:));
yfit = log(N_eps(:));

% Choose a fitting window to avoid:
% - too-large eps (coarse, not fractal)
% - too-small eps (noise / finite data saturation)
Np = size(D,1);

not_saturated = (N_eps < 0.2*Np);     % exclude where N is too close to Np
not_coarse    = (N_eps > 50);         % exclude where only a few boxes are hit

fit_range = not_saturated & not_coarse;

% If fit_range is empty (rare), loosen thresholds a bit:
if nnz(fit_range) < 5
    fit_range = (N_eps < 0.5*Np) & (N_eps > 20);
end

p = polyfit(xfit(fit_range), yfit(fit_range), 1);
D_box = p(1);


% ---- (4) Plot results ----
fprintf('N_boxes at largest eps: %d\n', N_eps(1));
fprintf('N_boxes at smallest eps: %d\n', N_eps(end));
fprintf('Number of points used: %d\n', size(D,1));

xfit = log(1 ./ eps_list(:));
yfit = log(N_eps(:));

local_slope = diff(yfit) ./ diff(xfit);

figure;
plot(eps_list(1:end-1), local_slope, 'o-'); grid on;
set(gca,'XScale','log');
xlabel('\epsilon'); ylabel('local slope');
title('Local slope vs \epsilon (pick a plateau)');
figure;
plot(xfit, yfit, 'o-','LineWidth',1.2); grid on; hold on;
plot(xfit(fit_range), polyval(p, xfit(fit_range)), 'LineWidth', 2);
xlabel('log(1/\epsilon)');
ylabel('log N(\epsilon)');
title(sprintf('Box-counting plot (Estimated D \\approx %.3f)', D_box));

fprintf('Estimated box-counting dimension D ≈ %.4f\n', D_box);


