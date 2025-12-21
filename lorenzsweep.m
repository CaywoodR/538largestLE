%% ============================================================
%  Lorenz system: autocorrelation + delay-coordinate visualization
% ============================================================

clear; clc; close all;

%% -------------------------------
% 1. Lorenz system definition
%% -------------------------------
function dXdt = lorenz_rhs(~, X)
    sigma = 10;
    beta  = 8/3;
    rho   = 28;

    x = X(1); y = X(2); z = X(3);
    dXdt = [ sigma*(y-x);
             x*(rho-z) - y;
             x*y - beta*z ];
end

%% -------------------------------
% 2. Time integration
%% -------------------------------
t0 = 0;
tf = 60;
dt = 0.01;                    % ← sampling interval (fixed!)
tspan = t0:dt:tf;

X0 = [0; 1; 1.05];
sol = ode45(@lorenz_rhs, [t0 tf], X0);

X = deval(sol, tspan).';       % sampled trajectory
t = tspan(:);

x = X(:,1);                   % scalar signal for embedding

%% -------------------------------
% 3. Remove transient (important!)
%% -------------------------------
burn_time = 10;               % seconds
burn_idx  = round(burn_time/dt);

x = x(burn_idx:end);
t = t(burn_idx:end);

x = x(:);
x = x - mean(x);              % demean for autocorr

dt_data = dt;                 % true sampling interval

%% -------------------------------
% 4. Autocorrelation
%% -------------------------------
maxLag = 300;
rho = zeros(maxLag+1,1);

den = sum(x.^2);
for k = 0:maxLag
    rho(k+1) = sum( x(1:end-k) .* x(1+k:end) ) / den;
end

lags = (0:maxLag).';
lags_time = lags*dt_data;

figure;
plot(lags_time, rho, 'k', 'LineWidth', 1.5, Color=[0 0.4470 0.7410]);
grid on;
xlabel('t (seconds)');
ylabel('A(\tau)');
title('Autocorrelation of Lorenz Attractor: X_{Time-series}');
yline(exp(-1),'--','1/e', Color= 'red');
yline(.123622,'--','Local Min', Color= 'blue');
yline(0,'-');

%% -------------------------------
% 5. Choose τ values in PHYSICAL TIME
%% -------------------------------
% Pick delays that increase monotonically in seconds
delay_times = [0.05 0.1 0.2 0.3 0.5 1.0];    % seconds
tauCandidates = round(delay_times / dt_data);

%% -------------------------------
% 6. 2D delay-coordinate visualization
%% -------------------------------
figure;
for i = 1:length(tauCandidates)
    tau = tauCandidates(i);

    if tau < length(x)
        subplot(2, ceil(length(tauCandidates)/2), i);
        plot(x(1:end-tau), x(1+tau:end), ...
             '.', 'MarkerSize', 2, 'Color', [0 0.4470 0.7410]);
        grid on;
        axis tight;

        title(sprintf('\\tau = %d  (%.2f s)', ...
              tau, tau*dt_data));
        xlabel('x(t)');
        ylabel('x(t+\tau)');
    end
end
sgtitle('2D delay plots for candidate delays \tau');

%% -------------------------------
% 7. Single selected τ example
%% -------------------------------
tau = tauCandidates(3);   % e.g. 0.2 s

figure;
plot(x(1:end-tau), x(1+tau:end), '.r', 'MarkerSize', 3);
grid on;
xlabel('x(t)');
ylabel(sprintf('x(t+\\tau), \\tau = %d samples', tau));
title(sprintf('2D delay plot at \\Delta t = %.2f s', tau*dt_data));
