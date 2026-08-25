function [feature_vec, feature_info] = ecg_amp_to_int_ratio_features(data, position, fs)

    % function [feature_vec, feature_info] = ecg_amp_to_int_ratio_features(data, position, fs)
    % Extract features related to the ECG amplitude-to- time interval ratios
    %
    % Inputs:
    %   data: ECG signal as a vector (in microvolts (uv)) (1D array).
    %   position: Fiducial points of the ECG signal (expressed as sample points).
    %   fs: Sampling frequency of the ECG signal (in Hz).
    %
    % Outputs:
    %   feature_vec contains features that represent ratios between ECG amplitudes and intervals.
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

    %% Define ECG points
    p = position.P;
    qrs_onset = position.QRSon;
    rpeak = position.R;
    qrs_offset = position.QRSoff;
    t = position.T;

    %% Calculate intervals
    pr_interval = convert_s_ms * (rpeak - p) / fs;
    qr_interval = convert_s_ms * (rpeak - qrs_onset) / fs;
    rs_interval = convert_s_ms * (qrs_offset - rpeak) / fs;
    rt_interval = convert_s_ms * ( t - rpeak) / fs;

    %% Calculate amplitude
    % PR amplitude
    pr_amp = NaN(size(rpeak)); 
    valid = ~isnan(rpeak) & ~isnan(p);
    pr_amp(valid) = data(rpeak(valid)) - data(p(valid));

    % QR amplitude
    qr_amp = NaN(size(rpeak)); 
    valid = ~isnan(rpeak) & ~isnan(qrs_onset);
    qr_amp(valid) = data(rpeak(valid)) - data(qrs_onset(valid));

    % RS amplitude
    rs_amp = NaN(size(rpeak)); 
    valid = ~isnan(rpeak) & ~isnan(qrs_offset);
    rs_amp(valid) = data(rpeak(valid)) - data(qrs_offset(valid));

    % RT amplitude
    rt_amp = NaN(size(rpeak)); 
    valid = ~isnan(rpeak) & ~isnan(t);
    rt_amp(valid) = data(rpeak(valid)) - data(t(valid));

    %% Calculate ratios
    ratio_pr = pr_amp./ pr_interval;
    ratio_qr = qr_amp./ qr_interval;
    ratio_rs = rs_amp./ rs_interval;
    ratio_rt = rt_amp./ rt_interval;

    % Omit NaN values
    valid_ratio_pr = ratio_pr(isfinite(ratio_pr));
    valid_ratio_qr = ratio_qr(isfinite(ratio_qr));
    valid_ratio_rs = ratio_rs(isfinite(ratio_rs));
    valid_ratio_rt = ratio_rt(isfinite(ratio_rt));

    %% Check if any valid values are present for valid_ratio_pr
    if isempty(valid_ratio_pr)
        mean_pr_ratios = NaN;
        std_pr_ratios = NaN;
        median_pr_ratios = NaN;
    else
        mean_pr_ratios = mean(valid_ratio_pr, "omitnan");
        std_pr_ratios = std(valid_ratio_pr, "omitnan");
        median_pr_ratios = median(valid_ratio_pr, "omitnan");
    end

    %% Check if any valid values are present for valid_ratio_qr
    if isempty(valid_ratio_qr)
        mean_qr_ratios = NaN;
        std_qr_ratios = NaN;
        median_qr_ratios = NaN;
    else
        mean_qr_ratios = mean(valid_ratio_qr, "omitnan");
        std_qr_ratios = std(valid_ratio_qr, "omitnan");
        median_qr_ratios = median(valid_ratio_qr, "omitnan");
    end
    
    %% Check if any valid values are present for valid_ratio_rs
    if isempty(valid_ratio_rs)
        mean_rs_ratios = NaN;
        std_rs_ratios = NaN;
        median_rs_ratios = NaN;
    else
        mean_rs_ratios = mean(valid_ratio_rs, "omitnan");
        std_rs_ratios = std(valid_ratio_rs, "omitnan");
        median_rs_ratios = median(valid_ratio_rs, "omitnan");
    end

    %% Check if any valid values are present for valid_ratio_rt
    if isempty(valid_ratio_rt)
        mean_rt_ratios = NaN;
        std_rt_ratios = NaN;
        median_rt_ratios = NaN;
    else
        mean_rt_ratios = mean(valid_ratio_rt, "omitnan");
        std_rt_ratios = std(valid_ratio_rt, "omitnan");
        median_rt_ratios = median(valid_ratio_rt, "omitnan");
    end

    %%  Results
    feature_vec = [mean_pr_ratios, std_pr_ratios, median_pr_ratios, ...
                            mean_qr_ratios, std_qr_ratios, median_qr_ratios,...
                            mean_rs_ratios, std_rs_ratios, median_rs_ratios,...
                            mean_rt_ratios, std_rt_ratios, median_rt_ratios];

    feature_vec = round(feature_vec,3);
    
    % Convert Inf to NaN
    feature_vec(isinf(feature_vec)) = NaN;

    % Define feature info
    feature_info.names = {"mean_pr_ratios", "std_pr_ratios", "median_pr_ratios", ...
                          "mean_qr_ratios", "std_qr_ratios", "median_qr_ratios", ...
                          "mean_rs_ratios", "std_rs_ratios", "median_rs_ratios", ...
                          "mean_rt_ratios", "std_rt_ratios", "median_rt_ratios"};

    feature_info.units = repmat({"uv/ms"}, 1, length(feature_info.names));
    feature_info.description = { ...
        "Mean of PR amplitude-to-interval ratios", ...
        "Standard deviation of PR amplitude-to-interval ratios", ...
        "Median of PR amplitude-to-interval ratios", ...
        "Mean of QR amplitude-to-interval ratios", ...
        "Standard deviation of QR amplitude-to-interval ratios", ...
        "Median of QR amplitude-to-interval ratios", ...
        "Mean of RS amplitude-to-interval ratios", ...
        "Standard deviation of RS amplitude-to-interval ratios", ...
        "Median of RS amplitude-to-interval ratios", ...
        "Mean of RT amplitude-to-interval ratios", ...
        "Standard deviation of RT amplitude-to-interval ratios", ...
        "Median of RT amplitude-to-interval ratios" ...
        };
end
