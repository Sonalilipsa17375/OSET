function [feature_vec, feature_info, mean_beat] = ecg_snr_features(data, R_peaks_indexes, increasing_beat_length, Pre_Rpeak_segment, Post_Rpeak_segment)
    
    % function [feature_vec, feature_info, mean_beat] = ecg_snr_features(data, R_peaks_indexes, 
    %                                                                                                             increasing_beat_length, Pre_Rpeak_segment,  Post_Rpeak_segment)
    % Extract features from ECG using beat signal-to-noise ratio (SNR)
    %
    % Inputs:
    %   data: ECG signal as a vector (in microvolts) (1D array).
    %   R_peaks_indexes:  A vector containing the R-peak indices of the ECG signal (expressed as sample points).
    %   increasing_beat_length
    %   Pre_Rpeak_segment: Portion used before the R-peak (float)
    %   Post_Rpeak_segment: Portion used after the R-peak (float)
    %
    % Outputs:
    % feature_vec: A vector contains
    %   1. median SNR
    %   2. mean SNR
    % feature_info: A structure contains feature descriptions, names and units.
    % mean_beat: Average ECG beat signal (1D array)
    %
    % Dependencies:
    %   1. `events_snr` function from the OSET package
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
    try 
            % Calculate beat length
            beat_length = round(increasing_beat_length * median(diff(R_peaks_indexes)));
                
            if mod(beat_length, 2) == 0
                 beat_length = beat_length + 1;
            end
                
            beat_length_modified = ceil([Pre_Rpeak_segment*beat_length, Post_Rpeak_segment*beat_length]);
                
            % Calculate SNR values
            [snr_median, snr_mean, mean_beat] = events_snr(data, R_peaks_indexes, beat_length_modified);
                
            % Assign the median and mean SNR values to the struct
            snr_median = snr_median(isfinite(snr_median));
            snr_mean = snr_mean(isfinite(snr_mean));
            snr_features.median = median(snr_median, 'omitnan');
            snr_features.mean = mean(snr_mean, 'omitnan');
                
            feature_vec = [snr_features.median, snr_features.mean];
            feature_vec = round(feature_vec, 3);
            % Convert Inf to NaN
            feature_vec(isinf(feature_vec)) = NaN;

    catch

            feature_vec = [NaN, NaN];
            mean_beat = nan(1, beat_length);

    end

    % Define feature info
    feature_info.names = {"snr_median", "snr_mean"};
    feature_info.units = {"dB", "dB"};
    feature_info.description = {"Median beat SNR", "Mean beat SNR"};

end
