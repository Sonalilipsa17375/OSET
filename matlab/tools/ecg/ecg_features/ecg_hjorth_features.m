function [feature_vec, feature_info] = ecg_hjorth_features(data, fs)

% [feature_vec, feature_info] = ecg_hjorth_features(data, fs)
% Extract features related to the complexity analysis of ECG signal
%
% Inputs:
%   data: ECG signal as a vector (in microvolts (uv)) (1D array).
%   fs:  Sampling frequency (Hz)
%
% Outputs:
%   feature_vec contains three features (activity, mobility and complexity)
%   feature_info contains feature descriptions, names and units
%
% Note: Derivatives are scaled by fs to match Hjorth's definition,
%       so Mobility is a true mean frequency in Hz (independent
%       of sampling rate). Per-sample implementations that omit fs report
%       Mobility in 1/sample and will differ from this by a factor of fs.
%
% Authors:
%   Seyedeh Somayyeh Mousavi
%   Sajjad Karimi
%   Reza Sameni
%
% bmemousavi@gmail.com
% Aug 2026
% Emory University, Georgia, USA

%%
% Calculate derivatives
dx = fs*diff(data);         % First derivative
ddx = fs*diff(dx);          % Second derivative

% Calculate variance
x_var = var(data);        % Activity
dx_var = var(dx);         % Variance of first derivative
ddx_var = var(ddx);     % Variance of second derivative

% Mobility and complexity calculations
if x_var == 0 || dx_var == 0
    mob = NaN;
    com = NaN;
else
    mob = sqrt(dx_var / x_var);                   % Mobility
    com = sqrt(ddx_var / dx_var) / mob;     % Complexity
end

feature_vec = [x_var, mob, com];
feature_vec = round(feature_vec ,3);
% Convert Inf to NaN
feature_vec(isinf(feature_vec)) = NaN;

% Define feature info
feature_info.names = {"Activity", "Mobility", "Complexity"};
feature_info.units = {"uv^2", "Hz", "scalar"};
feature_info.description = {"Hjorth parameters: Activity", ...
    "Hjorth parameters: Mobility", ...
    "Hjorth parameters: Complexity"};

end