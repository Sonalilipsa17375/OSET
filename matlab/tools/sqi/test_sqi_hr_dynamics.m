clear
close all
clc

data = importdata('/Users/rsameni/Downloads/FHR_values_per_recording.csv');
fs = 1000;
wlen = 10; % number of beats
hr = data(:, 2);
hr(hr <= 0) = [];
hr(~isfinite(hr)) = [];

rr = 60./ hr;
rr_samples = round(rr * fs);

peak_indexes = cumsum(rr_samples); 

[rho, mov_rho] = sqi_hr_dynamics(peak_indexes, fs, wlen);
mov_rho = abs(mov_rho);

mov_rho_norm = (mov_rho - min(mov_rho)) ./ (max(mov_rho) - min(mov_rho));

figure
subplot(211)
plot(rr)
grid

subplot(212)
plot(mov_rho)
grid

figure
% Create a vertical color-coded background
hold on
x = 1:length(mov_rho); % X values
yLimits = [min(rr), max(rr)]; % Y limits for the plot

% Loop to fill vertical strips for each mov_rho value
for i = 1:length(mov_rho)
    % Choose color for the current strip based on mov_rho
    % color = [0.7, 0.5, 0.5] * (mov_rho_norm(i)); % Adjust grayscale based on mov_rho
    color = [mov_rho_norm(i), 0.1, 1 - mov_rho_norm(i)]; % Adjust grayscale based on mov_rho
    fill([x(i), x(i), x(i)+1, x(i)+1], [yLimits(1), yLimits(2), yLimits(2), yLimits(1)], ...
         color, 'EdgeColor', 'none');
end

% Overlay the RR signal on top of the background
plot(rr, 'g', 'LineWidth', 1.5) % Plot RR in black
grid on
title('RR with Vertical Color-Coded Background (Signal Quality)')
xlabel('Index')
ylabel('RR')
xlim([1 length(mov_rho)]) % Set x-axis limits to match mov_rho
ylim(yLimits) % Set y-axis limits to match RR
% colorbar % Add a colorbar to explain mov_rho scale
% colormap(jet) % Use the jet colormap for color scaling

hold off