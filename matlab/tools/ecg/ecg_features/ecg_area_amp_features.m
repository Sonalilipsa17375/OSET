function [feature_vec, feature_info] = ecg_area_amp_features(data, position, fs)
        
    % function [feature_vec, feature_info] = ecg_area_amp_features(data, position, fs)
    % Extract features related to the ECG amplitude and area under the
    % curve of different ECG parts.
    %
    % Inputs:
    %   data: ECG signal as a vector (in microvolts (uv)) (1D array).
    %   position: Fiducial points of the ECG signal (expressed as sample points).
    %   fs: Sampling frequency of the ECG signal (in Hz).
    %
    % Outputs:
    %   feature_vec contains features related to: 
    %   1- ECG waves amplitudes  
    %   2- Area under the curve of the QRS complex, T wave and P wave,
    %   3- Amplitude of the ST segment, 
    %   4- Amplitude ratio of the QRS peak to the T waves.
    %    
    %   feature_info contains feature descriptions, names and units.
    %
    % Authors:
    %   Seyedeh Somayyeh Mousavi 
    %   Sajjad Karimi
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA

    %% Constant value
    convert_s_ms =1000;
    % Calculate the duration of each time step
    time_step = 1 / fs; 

    %% Define ECG points
    p_onset = position.Pon;
    p = position.P;
    p_offset = position.Poff;
    qrs_onset = position.QRSon;
    qrs = position.R;
    qrs_offset = position.QRSoff;
    t_onset = position.Ton;
    t = position.T;
    t_offset = position.Toff;
    N = length(data);

    %% Calculate amplitude
    % Extract T-peak amplitude from the ECG signal
    valid_indices_t = ~isnan(t) & ~isnan(t_onset) & ~isnan(t_offset) & ...
                                                        t >= 1 & t <= N & ...
                                                        t_onset >= 1 & t_onset <= N & ...
                                                        t_offset >= 1 & t_offset <= N;

    % Check if any valid amplitudes are present for T-wave
    if ~any(valid_indices_t)
        mean_t_amplitude = NaN;
        std_t_amplitude = NaN;
        median_t_amplitude = NaN;
        mean_t_auc = NaN;
        std_t_auc = NaN;
        median_t_auc = NaN;
        mean_t_abs_auc = NaN;
        std_t_abs_auc = NaN;
        median_t_abs_auc = NaN;
    
    else
        t_valid = t(valid_indices_t);
        t_onset_valid = t_onset(valid_indices_t);
        t_offset_valid = t_offset(valid_indices_t);   
    
        % Calculate statistics for T-wave
        amplitude_t_peak = data(t_valid) - (data(t_onset_valid) + data(t_offset_valid))/2;
        mean_t_amplitude = mean(amplitude_t_peak, "omitnan");
        std_t_amplitude = std(amplitude_t_peak, "omitnan");
        median_t_amplitude = median(amplitude_t_peak, "omitnan");
         
        % Compute AUC 
        t_auc = nan(size(t_valid)); 
        t_abs_auc = nan(size(t_valid));
        for i = 1:length(t_valid)
                
                t_auc(i) = sum(data(t_onset_valid(i):t_offset_valid(i))) * time_step * convert_s_ms;
                t_abs_auc(i) = sum(abs(data(t_onset_valid(i):t_offset_valid(i)))) * time_step * convert_s_ms;

        end
        mean_t_auc = mean(t_auc, "omitnan");
        std_t_auc = std(t_auc, "omitnan");
        median_t_auc = median(t_auc, "omitnan");
        mean_t_abs_auc = mean(t_abs_auc, "omitnan");
        std_t_abs_auc = std(t_abs_auc, "omitnan");
        median_t_abs_auc = median(t_abs_auc, "omitnan");
    end

    %% Calculate amplitude
    % Extract P-peak amplitude from the ECG signal
    valid_indices_p = ~isnan(p) & ~isnan(p_onset) & ~isnan(p_offset) & ...
                                                      p >= 1 & p <= N & ...
                                                      p_onset >= 1 & p_onset <= N & ...
                                                      p_offset >= 1 & p_offset <= N;

    % Check if any valid amplitudes are present for P-wave
    if ~any(valid_indices_p)
        mean_p_amplitude = NaN;
        std_p_amplitude = NaN;
        median_p_amplitude = NaN;
        mean_p_auc = NaN;
        std_p_auc = NaN;
        median_p_auc = NaN;
        mean_p_abs_auc = NaN;
        std_p_abs_auc = NaN;
        median_p_abs_auc = NaN;
    else
        % valid_indices_p = ~isnan(p) & ~isnan(p_onset) & ~isnan(p_offset);
        p_valid = p(valid_indices_p);
        p_onset_valid = p_onset(valid_indices_p);
        p_offset_valid = p_offset(valid_indices_p); 

        % Calculate statistics for P-wave
        amplitude_p_peak = data(p_valid) - (data(p_onset_valid) + data(p_offset_valid))/2;
        mean_p_amplitude = mean(amplitude_p_peak, "omitnan");
        std_p_amplitude = std(amplitude_p_peak, "omitnan");
        median_p_amplitude = median(amplitude_p_peak, "omitnan");

        % Compute AUC
        p_auc = nan(size(p_valid)); 
        p_abs_auc = nan(size(p_valid)); 

        for i = 1:length(p_valid)
                p_auc(i) = sum(data(p_onset_valid(i):p_offset_valid(i)))* time_step * convert_s_ms;
                p_abs_auc(i) = sum(abs(data(p_onset_valid(i):p_offset_valid(i))))* time_step * convert_s_ms;
        end
        
        mean_p_auc = mean(p_auc, "omitnan");
        std_p_auc = std(p_auc, "omitnan");
        median_p_auc = median(p_auc, "omitnan");
        mean_p_abs_auc = mean(p_abs_auc, "omitnan");
        std_p_abs_auc = std(p_abs_auc, "omitnan");
        median_p_abs_auc = median(p_abs_auc, "omitnan");

    end

    %% Extract QRS-peak amplitude from the ECG signal
    valid_indices_qrs = ~isnan(qrs) & ~isnan(qrs_onset) & ~isnan(qrs_offset) & ...
                                                        qrs >= 1 & qrs <= N & ...
                                                        qrs_onset >= 1 & qrs_onset <= N & ...
                                                        qrs_offset >= 1 & qrs_offset <= N;    
    
    % Check if any valid amplitudes are present for QRS complex
    if ~any(valid_indices_qrs)
        mean_qrs_amplitude = NaN;
        std_qrs_amplitude = NaN;
        median_qrs_amplitude = NaN;
        ratio_qrs_t_amplitude = NaN;
        ratio_qrs_p_amplitude = NaN;
        mean_qrs_auc = NaN;
        std_qrs_auc = NaN;
        median_qrs_auc = NaN;
        mean_qrs_abs_auc = NaN;
        std_qrs_abs_auc = NaN;
        median_qrs_abs_auc = NaN;
    else
        qrs_valid = qrs(valid_indices_qrs);
        qrs_onset_valid = qrs_onset(valid_indices_qrs);
        qrs_offset_valid = qrs_offset(valid_indices_qrs);

        % Calculate statistics for QRS complex
        amplitude_qrs_peak = data(qrs_valid) - (data(qrs_onset_valid) + data(qrs_offset_valid))/2;
        mean_qrs_amplitude = mean(amplitude_qrs_peak, "omitnan");
        std_qrs_amplitude = std(amplitude_qrs_peak, "omitnan");
        median_qrs_amplitude = median(amplitude_qrs_peak, "omitnan");
        ratio_qrs_t_amplitude = median_qrs_amplitude / median_t_amplitude;
        ratio_qrs_p_amplitude = median_qrs_amplitude / median_p_amplitude;

        % Compute AUC
        qrs_auc = nan(size(qrs_valid)); 
        qrs_abs_auc = nan(size(qrs_valid)); 

        for i = 1:length(qrs_valid)

                qrs_auc(i) = sum(data(qrs_onset_valid(i):qrs_offset_valid(i))) * time_step * convert_s_ms;
                qrs_abs_auc(i) = sum(abs(data(qrs_onset_valid(i):qrs_offset_valid(i)))) * time_step * convert_s_ms;

        end

        mean_qrs_auc = mean(qrs_auc, "omitnan");
        std_qrs_auc = std(qrs_auc, "omitnan");
        median_qrs_auc = median(qrs_auc, "omitnan");
        mean_qrs_abs_auc = mean(qrs_abs_auc, "omitnan");
        std_qrs_abs_auc = std(qrs_abs_auc, "omitnan");
        median_qrs_abs_auc= median(qrs_abs_auc, "omitnan");
    end    
    %% Calculate amplitude
    % Extract ST segment amplitude from the ECG signal
    valid_indices_st =  ~isnan(t_onset) & ~isnan(qrs_offset) & ...
                                                        t_onset >= 1 & t_onset <= N & ...
                                                        qrs_offset >= 1 & qrs_offset <= N ;   

    % Check if any valid amplitudes are present for the ST segment
    if ~any(valid_indices_st)
        mean_st_amplitude = NaN;
        std_st_amplitude = NaN;
        median_st_amplitude = NaN;    
    else
        st_onset_valid = t_onset(valid_indices_st);
        qrs_offset_valid = qrs_offset(valid_indices_st);
        amplitude_st = nan(size(st_onset_valid)); 
    
        % Extract amplitudes for the valid ST segment range
        for idx = 1:length(st_onset_valid)

            amplitude_st(idx) =  mean(data(qrs_offset_valid(idx):st_onset_valid(idx)), "omitnan");
            
        end

        % Calculate statistics 
        mean_st_amplitude = mean(amplitude_st, "omitnan");
        std_st_amplitude = std(amplitude_st, "omitnan");
        median_st_amplitude = median(amplitude_st, "omitnan");
    end
    
    %% Results
    feature_vec = [mean_qrs_amplitude, std_qrs_amplitude, median_qrs_amplitude, ...
                            mean_qrs_auc,           std_qrs_auc ,          median_qrs_auc, ...
                            mean_qrs_abs_auc,   std_qrs_abs_auc,   median_qrs_abs_auc,...
                            mean_t_amplitude,     std_t_amplitude,     median_t_amplitude, ...
                            ratio_qrs_t_amplitude, ...
                            mean_t_auc,               std_t_auc,               median_t_auc, ...
                            mean_t_abs_auc,       std_t_abs_auc,       median_t_abs_auc, ...
                            mean_p_amplitude,    std_p_amplitude,    median_p_amplitude, ...
                            ratio_qrs_p_amplitude, ...
                            mean_p_auc,              std_p_auc,              median_p_auc, ...
                            mean_p_abs_auc,      std_p_abs_auc,      median_p_abs_auc, ...
                            mean_st_amplitude,   std_st_amplitude,    median_st_amplitude];

    feature_vec = round(feature_vec ,3);
    % Convert Inf to NaN
    feature_vec(isinf(feature_vec)) = NaN;

    feature_info.names = {...
                "mean_qrs_amp", "std_qrs_amp", "median_qrs_amp", ...
                "mean_qrs_area", "std_qrs_area", "median_qrs_area", ...
                "mean_qrs_abs_area", "std_qrs_abs_area", "median_qrs_abs_area", ...
                "mean_t_amp", "std_t_amp", "median_t_amp", "ratio_qrs_t_amp", ...
                "mean_t_area", "std_t_area", "median_t_area",  ...
                "mean_t_abs_area", "std_t_abs_area", "median_t_abs_area", ...
                "mean_p_amp", "std_p_amp", "median_p_amp", "ratio_qrs_p_amp", ...
                "mean_p_area", "std_p_area", "median_p_area", ...
                "mean_p_abs_area", "std_p_abs_area", "median_p_abs_area", ...
                "mean_st_amp", "std_st_amp", "median_st_amp"
                };

    feature_info.units = { 
                "uv", "uv", "uv", ...
                "uv*ms", "uv*ms", "uv*ms", ...
                "uv*ms", "uv*ms", "uv*ms", ...
                "uv", "uv", "uv", "scalar", ...
                "uv*ms", "uv*ms", "uv*ms", ...
                "uv*ms", "uv*ms", "uv*ms",...
                "uv", "uv", "uv", "scalar", ...
                "uv*ms", "uv*ms", "uv*ms", ...
                "uv*ms", "uv*ms", "uv*ms",...
                "uv", "uv", "uv" };

    feature_info.description = { ...
                "Mean R-peak amplitude", ...
                "Standard deviation R-peak amplitude", ...
                "Median R-peak amplitude", ...
                "Mean QRS area", ...
                "Standard deviation QRS area", ...
                "Median QRS area", ...
                "Mean area of absolute QRS", ...
                "Standard deviation area of absolute QRS", ...
                "Median area of absolute QRS", ...
                "Mean T-peak amplitude", ...
                "Standard deviation T-peak amplitude", ...
                "Median T-peak amplitude", ...
                "QRS peak to T peak ratio", ...
                "Mean T wave area", ...
                "Standard deviation T wave area", ...
                "Median T wave area", ...
                "Mean area of absolute T wave", ...
                "Standard deviation area of absolute T wave", ...
                "Median area of absolute T wave", ...
                "Mean P-peak amplitude", ...
                "Standard deviation P-peak amplitude", ...
                "Median P-peak amplitude", ...
                "QRS peak to P peak ratio", ...
                "Mean P wave area", ...
                "Standard deviation P wave area", ...
                "Median P wave area", ...
                "Mean area of absolute P wave", ...
                "Standard deviation area of absolute P wave", ...
                "Median area of absolute P wave", ...
                "Mean ST-segment amplitude", ...
                "Standard deviation ST-segment amplitude", ...
                "Median ST-segment amplitude"
                 };

end
