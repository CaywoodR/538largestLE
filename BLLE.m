% estimate largest Lyapunov exponent from a time series


%data = readmatrix("Henon_Attractor.txt");

%data = readmatrix("Henon_delay_coords.txt");

data = readmatrix("Lorenz_delay_coords.txt");

%data = readmatrix("Lorenz_attractor_x.txt")';



% estimate largest Lyapunov exponent for the attractor
rng(42);
epsilon = 0.5;
max_delta_n = 20;








% choose some initial points randomly
s0 = round(rand(1,5).* (length(data)-max_delta_n));
theiler = 5;
sum_avgs = 0;

S_delta_n = zeros(max_delta_n, 1);

S = cell(length(s0), 1); % Initialize S as a cell array to store indices of close neighbors

for delta_n = 1:max_delta_n
    for ii = 1:length(s0)
        i0 = s0(ii);
        %disp(S{3});
        S{ii} = [];
        % build cell containing indices of close neighbors of current point
        for jj = 1:length(data)
            %if ((data(1,jj) - data(1,ii))^2 + (data(2,jj) - data(2,ii))^2)^(1/2) < epsilon && ii ~= jj
            if jj ~= i0 && abs(jj - i0) > theiler
            if norm(data(:,i0) - data(:, jj)) < epsilon
                S{ii} = [S{ii}, jj]; % Store indices in cell array
                
            end
            end
        end 
        % calculate the average distance from the current point after delta_n
        distances = zeros(length(S{ii}), 1);
        for kk = 1:length(S{ii})
            distances(kk) = norm(data(:, S{ii}(kk)+delta_n) - data(:, i0+delta_n)); % Compute distance to each neighbor
        end
        log_avg_dist = mean(log((distances)));
        sum_avgs = sum_avgs + log_avg_dist; % Accumulate the average log distances

    end
    S_delta_n(delta_n) = sum_avgs / length(s0); % Average log distance for this delta_n
    sum_avgs = 0; % Reset for the next delta_n

end

dt = 0.01;
xlow = 3;
xhi = 7;
t = (1:max_delta_n) * dt;

% plot S_delta_n
figure;
plot(t, S_delta_n, '-o'); 
hold on;

% plot fitted line over entire time axis
ys_fit = polyval(p, t);
plot(t, ys_fit, '--', 'LineWidth', 1);

xlabel('Time (t)');
ylabel('Average Log Distance');
title(['Average Log Distance vs time (t) -- Lyap. Exp. Estimate: ' num2str(slope)]);
grid on;



