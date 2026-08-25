function [feature_vec, feature_info] = ecg_time_intervals_features(position, fs)
    
    % function [feature_vec, feature_info] = ecg_time_intervals_features(position, fs)
    % Extract features related to the ECG time interval
    %
    % Inputs:
    %   position: Fiducial points of the ECG signal (expressed as sample points)
    %   fs: Sampling frequency (Hz)
    %
    % Outputs:
    % feature_vec contains features related to the ECG time intervals:
    %   - QRS complex duration
    %   - QT interval
    %   - PR interval
    %   - ST interval
    %   - PR segment
    %   - ST segment
    %   - Time intervals between peaks:
    %      * P to R peaks
    %      * Q to R peaks
    %      * S to R peaks
    %      * T to R peaks
    %   - Corrected QT intervals (QTC):
    %      * QTC_F (Fridericia correction)
    %      * QTC_B (Bazett correction)
    %
    % feature_info contains feature descriptions, names and units.
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
    convert_s_ms = 1000;
    convert_ms_s = 1/ convert_s_ms ;
    
    %% Define ECG points
    p_onset = position.Pon;
    p = position.P;
    p_offset = position.Poff;
    qrs_onset = position.QRSon;
    rpeak = position.R;
    qrs_offset = position.QRSoff;
    t_onset = position.Ton;
    t = position.T;
    t_offset = position.Toff;

    %% Calculate intervals
    rr_interval = convert_s_ms * diff(rpeak) / fs;
    mean_rr_interval = mean(rr_interval, "omitnan");

    qrs_complex = convert_s_ms * (qrs_offset - qrs_onset) / fs;
    qt_interval = convert_s_ms * (t_offset - qrs_onset) / fs;
    pr_interval = convert_s_ms * (qrs_onset - p_onset) / fs ;
    st_segment = convert_s_ms * (t_onset - qrs_offset) / fs ;    
    pr_segment = convert_s_ms * (qrs_onset - p_offset) / fs ;
    rt_peaks = convert_s_ms * ( t - rpeak) / fs;
    pr_peaks = convert_s_ms * (rpeak - p) / fs;
    qr_peaks = convert_s_ms * (rpeak - qrs_onset) / fs;
    rs_peaks = convert_s_ms * (qrs_offset - rpeak) / fs;

    % Omit NaN values
    valid_qrs_complex = qrs_complex(~isnan(qrs_complex) & qrs_complex >= 0);
    valid_qt_interval = qt_interval(~isnan(qt_interval) & qt_interval >= 0);
    valid_pr_interval = pr_interval(~isnan(pr_interval) & pr_interval >= 0);
    valid_st_segment = st_segment(~isnan(st_segment) & st_segment >= 0);
    valid_pr_segment = pr_segment(~isnan(pr_segment) & pr_segment >= 0);
    valid_rt_peaks = rt_peaks(~isnan(rt_peaks) & rt_peaks >= 0);
    valid_pr_peaks = pr_peaks(~isnan(pr_peaks) & pr_peaks >= 0);
    valid_qr_peaks = qr_peaks(~isnan(qr_peaks) & qr_peaks >= 0);
    valid_rs_peaks = rs_peaks(~isnan(rs_peaks) & rs_peaks >= 0);

    %% Check if any valid intervals are present for qrs_complex
    if isempty(valid_qrs_complex)
        mean_qrs = NaN;
        std_qrs = NaN;
        median_qrs = NaN;
    else
        % Calculate statistics
        mean_qrs = mean(valid_qrs_complex, "omitnan");
        std_qrs = std(valid_qrs_complex, "omitnan");
        median_qrs = median(valid_qrs_complex, "omitnan");
    end

    %% Check if any valid intervals are present for qt_interval
    if isempty(valid_qt_interval)
        mean_qt_interval = NaN;
        std_qt_interval = NaN;
        median_qt_interval = NaN;
    else
        % Calculate statistics
        mean_qt_interval = mean(valid_qt_interval, "omitnan");
        std_qt_interval = std(valid_qt_interval, "omitnan");
        median_qt_interval = median(valid_qt_interval, "omitnan");
    end

    %% QTc_b and QTc_f    
    if ~isnan(mean_rr_interval) && mean_rr_interval > 0 && ~isnan(mean_qt_interval)
            QTc_b = mean_qt_interval / ((mean_rr_interval*convert_ms_s)^(1/2));
            QTc_f = mean_qt_interval / ((mean_rr_interval*convert_ms_s)^(1/3));
    else
            QTc_b = NaN;
            QTc_f = NaN;
    end
    
    %% Check if any valid intervals are present for pr_interval
    if isempty(valid_pr_interval)
        mean_pr_interval = NaN;
        std_pr_interval = NaN;
        median_pr_interval = NaN;
    else
        % Calculate statistics for pr_interval
        mean_pr_interval = mean(valid_pr_interval, "omitnan");
        std_pr_interval = std(valid_pr_interval, "omitnan");
        median_pr_interval = median(valid_pr_interval, "omitnan");
    end

    %% Check if any valid intervals are present for st_segment
    if isempty(valid_st_segment)
        mean_st_segment= NaN;
        std_st_segment = NaN;
        median_st_segment = NaN;
    else
    
        % Calculate statistics for st_interval
        mean_st_segment = mean(valid_st_segment, "omitnan");
        std_st_segment = std(valid_st_segment, "omitnan");
        median_st_segment = median(valid_st_segment, "omitnan");
    end
    
    %% Check if any valid intervals are present for pr_segment
    if isempty(valid_pr_segment)
        mean_pr_segment = NaN;
        std_pr_segment = NaN;
        median_pr_segment = NaN;
    else
    
        % Calculate statistics for pr_segment
        mean_pr_segment = mean(valid_pr_segment, "omitnan");
        std_pr_segment = std(valid_pr_segment, "omitnan");
        median_pr_segment = median(valid_pr_segment, "omitnan");
    end

    %% Check if any valid intervals are present for valid_rt_peaks
    if isempty(valid_rt_peaks)
        mean_rt_peaks = NaN;
        std_rt_peaks = NaN;
        median_rt_peaks = NaN;
    else
    
        % Calculate statistics for valid_rt_peaks
        mean_rt_peaks = mean(valid_rt_peaks, "omitnan");
        std_rt_peaks = std(valid_rt_peaks, "omitnan");
        median_rt_peaks = median(valid_rt_peaks, "omitnan");
    end

    %% Check if any valid intervals are present for valid_pr_peaks
    if isempty(valid_pr_peaks)
        mean_pr_peaks = NaN;
        std_pr_peaks = NaN;
        median_pr_peaks = NaN;
    else
    
        % Calculate statistics for valid_pr_peaks
        mean_pr_peaks = mean(valid_pr_peaks, "omitnan");
        std_pr_peaks = std(valid_pr_peaks, "omitnan");
        median_pr_peaks = median(valid_pr_peaks, "omitnan");
    end

    %% Check if any valid intervals are present for valid_qr_peaks
    if isempty(valid_qr_peaks)
        mean_qr_peaks = NaN;
        std_qr_peaks = NaN;
        median_qr_peaks = NaN;
    else
    
        % Calculate statistics for valid_qr_peaks
        mean_qr_peaks = mean(valid_qr_peaks, "omitnan");
        std_qr_peaks = std(valid_qr_peaks, "omitnan");
        median_qr_peaks = median(valid_qr_peaks, "omitnan");
    end
    
    %% Check if any valid intervals are present for valid_rs_peaks
    if isempty(valid_rs_peaks)
        mean_rs_peaks = NaN;
        std_rs_peaks = NaN;
        median_rs_peaks = NaN;
    else
    
        % Calculate statistics for valid rs_peaks
        mean_rs_peaks = mean(valid_rs_peaks, "omitnan");
        std_rs_peaks = std(valid_rs_peaks, "omitnan");
        median_rs_peaks = median(valid_rs_peaks, "omitnan");
    end

    feature_vec = [mean_pr_interval, std_pr_interval, median_pr_interval, ...
                           mean_qt_interval, std_qt_interval, median_qt_interval, ...
                           QTc_b, QTc_f, ...
                           mean_st_segment, std_st_segment, median_st_segment, ...
                           mean_pr_segment, std_pr_segment, median_pr_segment, ...
                           mean_pr_peaks, std_pr_peaks, median_pr_peaks, ...
                           mean_qr_peaks, std_qr_peaks, median_qr_peaks, ...
                           mean_qrs, std_qrs, median_qrs, ...
                           mean_rs_peaks, std_rs_peaks, median_rs_peaks, ...
                           mean_rt_peaks, std_rt_peaks, median_rt_peaks];

    % Convert Inf to NaN
    feature_vec(isinf(feature_vec)) = NaN;

    feature_vec = round(feature_vec, 3);
    
    % Define feature info
    feature_info.names = {"mean_pr_interval", "std_pr_interval", "median_pr_interval", ...
                                        "mean_qt_interval", "std_qt_interval", "median_qt_interval", ...
                                        "qtc_b", "qtc_f", ...
                                        "mean_st_segment", "std_st_segment", "median_st_segment", ...
                                        "mean_pr_segment", "std_pr_segment", "median_pr_segment", ...
                                        "mean_pr_peaks_interval", "std_pr_peaks_interval", "median_pr_peaks_interval", ...
                                        "mean_qr_peaks_interval", "std_qr_peaks_interval", "median_qr_peaks_interval", ...
                                        "mean_qrs_complex_interval", "std_qrs_complex_interval", "median_qrs_complex_interval", ...
                                        "mean_rs_peaks_interval", "std_rs_peaks_interval", "median_rs_peaks_interval", ...
                                        "mean_rt_peaks_interval", "std_rt_peaks_interval", "median_rt_peaks_interval"};

    feature_info.units = repmat({"ms"}, 1, length(feature_info.names));

    feature_info.description = {"Mean PR-interval", "Standard deviation of PR-interval", "Median PR-interval", ...
                                               "Mean QT-interval", "Standard deviation of QT-interval", "Median QT-interval", ...
                                               "QTc (Bazett)","QTc (Fridericia)", ...
                                               "Mean ST-segment (QRSoff-Ton)", "Standard deviation of ST-segment (QRSoff-Ton)", "Median ST-segment (QRSoff-Ton)", ...
                                               "Mean PR-segment (Poff-QRSon)", "Standard deviation of PR-segment (Poff-QRSon)", "Median PR-segment (Poff-QRSon)", ...
                                               "Mean P to R peaks interval", "Standard deviation of P to R peaks interval", "Median P to R peaks interval", ...
                                               "Mean QRSon to R peak interval", "Standard deviation of QRSon to R peak interval", "Median QRSon to R peak interval", ...
                                               "Mean QRS-complex interval", "Standard deviation of QRS-complex interval", "Median  QRS-complex interval", ...
                                               "Mean R peak to QRSoff interval", "Standard deviation of R peak to QRSoff interval", "Median R peak to QRSoff interval", ...
                                               "Mean R peak to T peak interval", "Standard deviation of R peak to T peak interval", "Median R peak to T peak interval"};
end
