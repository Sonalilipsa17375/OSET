function [feature_vec, feature_info] = ecg_svd_features(data, R_peak_indexes, number_eigenvalues)
    
    % function [feature_vec, feature_info] = ecg_svd_features(data, R_peak_indexes, number_eigenvalues)
    % Extract features from ECG using Singular Value Decomposition (SVD)
    %
    % Inputs:
    %   data: ECG signal as a vector (in microvolts) (1D array).
    %   R_peaks_indexes - A vector containing the R-peak indices of the ECG signal (expressed as sample points).
    %   number_eigenvalues: Number of eigenvalues to consider for the feature vector (scaler).
    %
    % Outputs:
    %   feature_vec contains normalized singular values in percentage.
    %   feature_info contains feature descriptions, names and units.
    %
    %
    % Dependencies:
    %   1. `event_stacker` function from the OSET package
    %
    % Authors:
    %   Seyedeh Somayyeh Mousavi 
    %   Sajjad Karimi
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA

    %% Check
    if isempty(number_eigenvalues) || number_eigenvalues < 1
        number_eigenvalues = 1;
        warning('ecg_svd_features:invalidNumEigenvalues', ...
            'number_eigenvalues must be >= 1. Defaulting to 1.');
    end
    number_eigenvalues = round(number_eigenvalues);

    try
            % Step 1: Compute RR intervals in samples
            RR_intervals_samples = diff(R_peak_indexes);

            if numel(R_peak_indexes) < 2 || isempty(RR_intervals_samples) || all(isnan(RR_intervals_samples))
                        error('ecg_svd_features:insufficientPeaks', ...
                                'Need at least 2 valid R-peaks to compute RR intervals.');
            end
        
            % Step 2: Determine the event bounds (median of RR intervals)
            event_bounds = round(median(RR_intervals_samples));
        
            % Ensure event_bounds is odd
            if mod(event_bounds, 2) == 0

                    event_bounds = event_bounds + 1;

            end

        
            % Step 3: Extract signal segments using event_stacker
            [stacked_events, ~] = event_stacker(data, R_peak_indexes, event_bounds);
        
            % Step 4: Perform Singular Value Decomposition (SVD)
            singular_values = svd(stacked_events);
        
            % Step 5: Normalize singular values
            singular_values = singular_values / sum(singular_values);
        
            % Store the normalized singular values into the output array
            if length(singular_values) >= number_eigenvalues
                % Truncate if there are more eigenvalues than desired
                singularvalues_normalized = singular_values(1:number_eigenvalues);
            else
                % Zero-pad if there are fewer eigenvalues than desired
                singularvalues_normalized = [singular_values; zeros(number_eigenvalues - length(singular_values), 1)];
            end
        
            % Step 6: Convert normalized singular values to percentages
            feature_vec = singularvalues_normalized' * 100;
        
            feature_vec = round(feature_vec, 3);
            % Convert Inf to NaN
            feature_vec(isinf(feature_vec)) = NaN;

    catch

            feature_vec = nan(1, number_eigenvalues);
            
    end

    % Define feature info
    for n = 1:number_eigenvalues
        feature_info.names{n} = string(['svd', num2str(n)]);
        feature_info.units{n} = "scalar";
        feature_info.description{n} = string(['Normalized SVD value ', num2str(n)]);
    end

end

