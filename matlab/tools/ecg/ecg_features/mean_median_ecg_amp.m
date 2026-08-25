function [mean_amplitude, median_amplitude] = mean_median_ecg_amp(data, position, length_data)

    % function [mean_amplitude, median_amplitude] = mean_median_ecg_amp(data, position, length_data)
    % Extract features related to the ECG amplitude
    %
    % Inputs:
    %   data: ECG signal as a vector (in microvolts (uv)) (1D array).
    %   position: Fiducial points of the ECG signal (expressed as sample points).
    %   length_data
    % Outputs:
    %   mean_amplitude: mean ECG waves amplitudes 
    %   median_amplitude: median ECG waves amplitudes 
    %
    % Authors:
    %   Seyedeh Somayyeh Mousavi 
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA
    %%
    try
            %% Define ECG points and valid ones
            N = length_data;
            p_onset = position.Pon;
            p = position.P;
            p_offset = position.Poff;
        
            qrs_onset = position.QRSon;
            qrs = position.R;
            qrs_offset = position.QRSoff;
            
            t_onset = position.Ton;
            t = position.T;
            t_offset = position.Toff;
        
            valid_indices_p = ~isnan(p) & ~isnan(p_onset) & ~isnan(p_offset) & ...
                                                                p >= 1 & p <= N & ...
                                                                p_onset >= 1 & p_onset <= N & ...
                                                                p_offset >= 1 & p_offset <= N;
        
            valid_indices_t = ~isnan(t) & ~isnan(t_onset) & ~isnan(t_offset) & ...
                                                                t >= 1 & t <= N & ...
                                                                t_onset >= 1 & t_onset <= N & ...
                                                                t_offset >= 1 & t_offset <= N;
        
            valid_indices_qrs = ~isnan(qrs) & ~isnan(qrs_onset) & ~isnan(qrs_offset) & ...
                                                                qrs >= 1 & qrs <= N & ...
                                                                qrs_onset >= 1 & qrs_onset <= N & ...
                                                                qrs_offset >= 1 & qrs_offset <= N;   
        
            %% Calculate amplitude
            % Extract P-peak amplitude from the ECG signal
            if ~any(valid_indices_p)
        
                 mean_amplitude.P = NaN;
                 median_amplitude.P = NaN;
            
            else
                 p_valid = p(valid_indices_p);
                 p_onset_valid = p_onset(valid_indices_p);
                 p_offset_valid = p_offset(valid_indices_p);   
            
                 % Calculate statistics for P-wave
                 amplitude_p_peak = data(p_valid) - (data(p_onset_valid) + data(p_offset_valid))/2;
                 mean_amplitude.P = mean(amplitude_p_peak, "omitnan");
                 median_amplitude.P = median(amplitude_p_peak, "omitnan");
                 
            end
        
            % Extract QRS-peak amplitude from the ECG signal
            if ~any(valid_indices_qrs)
        
                 mean_amplitude.R = NaN;
                 median_amplitude.R = NaN;
            
            else
                 qrs_valid = qrs(valid_indices_qrs);
                 qrs_onset_valid = qrs_onset(valid_indices_qrs);
                 qrs_offset_valid = qrs_offset(valid_indices_qrs);   
            
                 % Calculate statistics for QRS-wave
                 amplitude_qrs_peak = data(qrs_valid) - (data(qrs_onset_valid) + data(qrs_offset_valid))/2;
                 mean_amplitude.R = mean(amplitude_qrs_peak, "omitnan");
                 median_amplitude.R = median(amplitude_qrs_peak, "omitnan");
                 
            end
        
            % Extract T-peak amplitude from the ECG signal
            if ~any(valid_indices_t)
        
                 mean_amplitude.T = NaN;
                 median_amplitude.T = NaN;
            
            else
                 t_valid = t(valid_indices_t);
                 t_onset_valid = t_onset(valid_indices_t);
                 t_offset_valid = t_offset(valid_indices_t);   
            
                 % Calculate statistics for T-wave
                 amplitude_t_peak = data(t_valid) - (data(t_onset_valid) + data(t_offset_valid))/2;
                 mean_amplitude.T = mean(amplitude_t_peak, "omitnan");
                 median_amplitude.T = median(amplitude_t_peak, "omitnan");
                 
            end
    catch

            mean_amplitude.P = NaN;   median_amplitude.P = NaN;
            mean_amplitude.R = NaN;   median_amplitude.R = NaN;
            mean_amplitude.T = NaN;   median_amplitude.T = NaN;
            
    end 
end