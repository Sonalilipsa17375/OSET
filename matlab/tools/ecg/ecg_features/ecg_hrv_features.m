function [feature_vec, feature_info] = ecg_hrv_features( data, R_peaks_indexes, fs, ...
                                                                            MinRPeaks, MinRRFactor, MaxRRFactor, ...
                                                                            MinRRSeconds, MaxRRSeconds, ...
                                                                            MinRRPercentile, MaxRRPercentile, ...
                                                                            MinHRPercentile, MaxHRPercentile)
    % function [feature_vec, feature_info] = ecg_hrv_features( data, R_peaks_indexes, fs, ...
    %                                                                   MinRPeaks, MinRRFactor, MaxRRFactor, ...
    %                                                                   MinRRSeconds, MaxRRSeconds, ...
    %                                                                   MinRRPercentile, MaxRRPercentile, ...
    %                                                                   MinHRPercentile, MaxHRPercentile, ...
    %                                                                   )
    % Extract Heart Rate Variability (HRV) features.
    %
    % Inputs:
    %   data: ECG signal in microvolts (1D array)
    %   R_peaks_indexes:  A vector containing the R-peak indices of the ECG signal 
    %   (expressed as sample points)
    %   fs:  Sampling frequency (Hz)
    %   MinRPeaks: Minimum number of R-peaks required (Default: 4)
    %   MinRRFactor: Minimum RR as a fraction of median RR (Default: 0.5)
    %   MaxRRFactor: Maximum RR as a fraction of median RR (Default: 2)
    %   MinRRSeconds: Absolute minimum RR interval (s) (Default: 0.2)
    %   MaxRRSeconds: Absolute maximum RR interval (s) (Default: 2.0)
    %   MinRRPercentile: Minimum range for RR filtering (Default: 2.5)
    %   MaxRRPercentile: Maximum range for RR filtering (Default: 97.5)
    %   MinHRPercentile: Minimum range for HR filtering (Default: 2.5)
    %   MaxHRPercentile: Maximum range for HR filtering (Default: 97.5)
    %
    % Outputs:
    % feature_vec contains the following HRV features (12):
    %   1. ECG_Length in seconds
    %   2. Number of beats
    %   3. RMSSD (Root Mean Square of Successive Differences) in milliseconds
    %   4. SDNN (Standard Deviation of normal-to-normal (NN) intervals) in milliseconds
    %   5. HR (MEDIAN, MEAN, upper_5, lower_5)
    %   6. HRF (hrf_pip, hrf_ials, hrf_pnn_ss, hrf_pnn_as)
    %
    % feature_info contains feature descriptions, names and units.
    %    
    % Dependencies:
    %   1. `heart_rate_fragmentation` function from the OSET package
    %
    % Authors:
    %   Seyedeh Somayyeh Mousavi 
    %   Sajjad Karimi
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA

%% Default values for optional parameters
if nargin < 4 || isempty(MinRPeaks)
    MinRPeaks = 4;
end

if nargin < 5 || isempty(MinRRFactor)
    MinRRFactor = 0.5;
end

if nargin < 6 || isempty(MaxRRFactor)
    MaxRRFactor = 2;
end

% Corresponds to HR = 60 / 0.2 = 300 bpm (upper bound, shortest allowed RR)
if nargin < 7 || isempty(MinRRSeconds)
    MinRRSeconds = 0.2;  
end

% Corresponds to HR = 60 / 2.0 = 30 bpm (lower bound, longest allowed RR)
if nargin < 8 || isempty(MaxRRSeconds)
    MaxRRSeconds = 2;
end

if nargin < 9 || isempty(MinRRPercentile)
    MinRRPercentile = 2.5;
end

if nargin < 10 || isempty(MaxRRPercentile)
    MaxRRPercentile = 97.5;
end

if nargin < 11 || isempty(MinHRPercentile)
    MinHRPercentile = 2.5;
end

if nargin < 12 || isempty(MaxHRPercentile)
    MaxHRPercentile = 97.5;
end

% Constant value
convert_s_ms =1000;

% Initialize with NaNs so any early failure still returns a well-formed struct
features.Length_data = length(data)/fs;
features.N_beats = length(R_peaks_indexes);
features.RMSSD = nan;
features.SDNN = nan;
features.HR_median = nan;
features.HR_mean = nan;
features.HR_upper_5  = nan;
features.HR_lower_5   = nan;
features.hrf_pip = nan;
features.hrf_ials  = nan;
features.hrf_pnn_ss  = nan;
features.hrf_pnn_as  = nan;

% Extract features if length(R_peaks) >= MinRPeaks
if length(R_peaks_indexes) >= MinRPeaks

    try     
            % Compute RR intervals in samples
            RR_intervals_samples = diff(R_peaks_indexes); % Time between successive R-peaks in samples
        
            % Convert RR intervals to time in seconds using sampling frequency (fs)
            RR_intervals_seconds = RR_intervals_samples / fs;
            
            % Filtering for non-normal beats
            med = median(RR_intervals_seconds);
            min_rr = max(MinRRFactor*med, MinRRSeconds);
            max_rr = min(MaxRRFactor*med, MaxRRSeconds);
            RR_intervals_seconds(RR_intervals_seconds < min_rr) = [];
            RR_intervals_seconds(RR_intervals_seconds > max_rr) = [];
        
            % Additional percentile-based outlier filtering.
            RR_percentiles = prctile(RR_intervals_seconds, [MinRRPercentile  MaxRRPercentile]);
            filtered_RR_intervals_seconds = RR_intervals_seconds(RR_intervals_seconds >= RR_percentiles(1) & ...
                                                                                                      RR_intervals_seconds <= RR_percentiles(2));
        
            % Feature 1: RMSSD (Root Mean Square of Successive Differences)
            dRR = diff(filtered_RR_intervals_seconds);
            % Convert to milliseconds    
            RMSSD = sqrt(mean(dRR.^2))* convert_s_ms;
        
            % Feature 2: SDNN (Standard Deviation of NN intervals)
            SDNN = std(filtered_RR_intervals_seconds);
            % Convert to milliseconds
            SDNN = SDNN * convert_s_ms;
        
            % Calculate heart rate (HR)
            HR = 60 ./ filtered_RR_intervals_seconds;  % HR in beats per minute based on RR intervals in seconds
            HR_percentiles = prctile(HR, [MinHRPercentile  MaxHRPercentile]);
            HR = HR(HR >= HR_percentiles(1) & HR <= HR_percentiles(2));  % Filter HR values within the range
            
            % Calculate mean heart rate
            mean_HR = mean(HR);
        
            % Calculate median heart rate
            median_HR = median(HR);
        
            % Calculate heart rate for lower 5% 
            HR_lower_5 = prctile(HR, 5);
        
            % Calculate heart rate for upper 5% 
            HR_upper_5 = prctile(HR, 95);
        
            % Calculate heart rate for the interquartile range (5%-95%)
            [pip, ials, pnn_ss, pnn_as ] = heart_rate_fragmentation(filtered_RR_intervals_seconds * convert_s_ms , fs);

            % Store the HRV features
            features.Length_data = length(data)/fs;
            features.N_beats = length(filtered_RR_intervals_seconds) + 1;
            features.RMSSD = RMSSD;
            features.SDNN = SDNN;
            features.HR_median = median_HR;
            features.HR_mean = mean_HR;
            features.HR_upper_5 = HR_upper_5;
            features.HR_lower_5 = HR_lower_5;
            features.hrf_pip = pip;
            features.hrf_ials = ials;
            features.hrf_pnn_ss = pnn_ss;
            features.hrf_pnn_as = pnn_as;

    catch ME
                warning('ecg_hrv_features:computationFailed', '%s', ME.message);
    end

end


feature_vec = [ ...
    round(features.Length_data, 3), ...
    features.N_beats, ...
    round(features.RMSSD, 3), ...
    round(features.SDNN, 3), ...
    round(features.HR_median, 2), ...
    round(features.HR_mean, 2), ...
    round(features.HR_upper_5, 2), ...
    round(features.HR_lower_5, 2), ...
    round(features.hrf_pip, 3), ...
    round(features.hrf_ials, 3), ...
    round(features.hrf_pnn_ss, 3), ...
    round(features.hrf_pnn_as, 3) ...
];

% Convert Inf to NaN
feature_vec(isinf(feature_vec)) = NaN;

% Define feature info
feature_info.names = {"ECG_len", "n_beats", ...
                                    "rmssd", "sdnn", ...
                                    "hr_median", "hr_mean", ...
                                    "hr_upper_5", "hr_lower_5", ...
                                    "hrf_pip", "hrf_ials", "hrf_pnn_ss", "hrf_pnn_as"};
feature_info.units =   {"s", "scalar", ...
                                   "ms", "ms", ...
                                   "bpm", "bpm", ...
                                   "bpm", "bpm", ...
                                   "scalar","scalar","scalar","scalar"};
feature_info.description = {"ECG length", "Number of ECG beats", ...
                                           "Root Mean Square of Successive Differences for HRV", ...
                                           "Standard Deviation of Normal-to-Normal Intervals", ...
                                           "Median heart-rate", "Mean heart-rate", ...
                                           "95% heart-rate", "5% heart-rate",...
                                           "Heart rate fragmentation: Percentage of inflection points", ...
                                           "Heart rate fragmentation: " + ...
                                           "Inverse average length of the acceleration/deceleration segments",...
                                           "Heart rate fragmentation: " + ...
                                           "Percentage of short segments", ...
                                           "Heart rate fragmentation: " + ...
                                           "Percentage of NN intervals in alternation segments"};

end
