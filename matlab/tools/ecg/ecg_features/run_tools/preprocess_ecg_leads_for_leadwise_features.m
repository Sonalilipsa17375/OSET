function data_filtered = preprocess_ecg_leads_for_leadwise_features(data, fs, f_notch, Q_factor, low_pass_fre, high_pass_fre, median_window_coff, mean_window_coff)

    % data_filtered = preprocess_ecg_leads_for_leadwise_features(data, fs, f_notch, 
    % Q_factor, low_pass_fre, high_pass_fre, median_window_coff, mean_window_coff)
    % Preprocess ECG signal: NaN handling, notch filtering, baseline removal, and
    % zero-phase low-pass/high-pass filtering.
    %
    % Inputs:
    %   data          (channels x samples double): raw ECG signal
    %   fs            (scalar): sampling frequency (Hz)
    %   f_notch       (scalar): notch filter frequency (Hz)
    %   Q_factor      (scalar): quality factor for the notch filter
    %   low_pass_fre  (scalar): low-pass cut-off frequency (Hz)
    %   high_pass_fre (scalar): high-pass cut-off frequency (Hz)
    %   median_window_coeff (scalar): coefficient for the moving median window
    %   mean_window_coeff (scalar): coefficient for the moving mean window
    %
    % Outputs:
    %   data_filtered: preprocessed signal same size as data.
    %
    % Author: 
    %   Seyedeh Somayyeh Mousavi
    %   Reza Sameni
    %   Emory University, Georgia, USA
    %   bmemousavi@gmail.com
    %
    % AUG, 2026
    
    %%
    n_channels = size(data, 1);

    % Design notch filter
    W0 = f_notch / (fs/2);
    [b, a] = iirnotch(W0, W0 / Q_factor);

    data_filtered = zeros(size(data));

    % Preprocessing
    for c = 1:n_channels
        ecg_raw = data(c,:);

        % Handle nan values
        ecg_raw(isnan(ecg_raw)) = median(ecg_raw,'omitnan');

        % Apply notch filter
        ecg_raw1 = filtfilt(b, a, ecg_raw);

        % Remove baseline
        ecg_raw2 = ecg_raw1 - (movmean(movmedian(ecg_raw1(:),[round(median_window_coff*fs),round(median_window_coff*fs)]),[round(mean_window_coff*fs),round(mean_window_coff*fs)]))';

        % Apply high-pass and low-pass filters (Zero-phase filtering)
        ecg_raw3 = ecg_raw2 - lp_filter_zero_phase(ecg_raw2, low_pass_fre / fs);
        ecg_raw4 = lp_filter_zero_phase(ecg_raw3, high_pass_fre / fs);

        data_filtered(c,:) = ecg_raw4;
    end
end