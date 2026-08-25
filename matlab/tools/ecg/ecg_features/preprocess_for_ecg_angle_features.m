function data_filtered = preprocess_for_ecg_angle_features(data, fs, low_pass_fre, high_pass_fre, median_window_coff, mean_window_coff)

    % data_filtered = preprocess_for_ecg_angle_features(data, fs, low_pass_fre, high_pass_fre, median_window_coff, mean_window_coff)
    % Preprocess ECG signal: NaN handling, notch filtering, baseline removal, and
    % zero-phase low-pass/high-pass filtering.
    %
    % Inputs:
    %   data          (channels x samples double): raw ECG signal
    %   fs            (scalar): sampling frequency (Hz)
    %   low_pass_fre  (scalar): low-pass cut-off frequency (Hz)
    %   high_pass_fre (scalar): high-pass cut-off frequency (Hz)
    %   median_window_coeff (scalar): coefficient for the moving median window
    %   mean_window_coeff (scalar): coefficient for the moving mean window
    %
    % Outputs:
    %   data_filtered: preprocessed signal same size as data.
    %
    % Authors: 
    %   Seyedeh Somayyeh Mousavi
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA
    
    %%
    data_filtered = zeros(size(data));
    n_channels = size(data, 1);
    % Preprocessing
    for c = 1:n_channels

                    ecg_raw = data(c,:);
                    
                    % Handle nan values 
                    ecg_raw(isnan(ecg_raw)) = median(ecg_raw,'omitnan');

                    %  Remove baseline 
                    ecg_raw1 = ecg_raw - (movmean(movmedian(ecg_raw(:),[round(median_window_coff*fs),round(median_window_coff*fs)]),[round(mean_window_coff*fs),round(mean_window_coff*fs)]))';

                    % Apply high-pass and low-pass filters (Zero-phase filtering)
                    ecg_raw2 = ecg_raw1 - lp_filter_zero_phase(ecg_raw1, low_pass_fre / fs);
                    ecg_raw3 = lp_filter_zero_phase(ecg_raw2, high_pass_fre / fs);

                    data_filtered(c,:) = ecg_raw3;
    end
        
end