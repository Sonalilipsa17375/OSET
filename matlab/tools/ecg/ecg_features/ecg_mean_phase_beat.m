function [feature_vec, feature_info] = ecg_mean_phase_beat(data, r_peak_impulse_train, num_phase_beans, Flag_baseline_alignment)

    % function [feature_vec, feature_info] = ecg_mean_phase_beat(data, r_peak_impulse_train, num_phase_beans, Flag_baseline_alignment)
    % Extract features related to the morphology of the average ECG
    % beat using phase domain information
    %
    % Inputs:
    %   data:  ECG signal as a vector (in microvolts (uV)) (1D array).
    %   r_peak_impulse_train: An impulse train vector containing 1s at R-peak locations of the ECG signal 
    %   (expressed as sample points).
    %   num_phase_beans: number of phase bins
    %   Flag_baseline_alignment (default = 1)
    %
    % Outputs:
    %   feature_vec contains the mean beat
    %   feature_info contains feature descriptions, names and units
    %
    % Dependencies:
    %   1- phase_calculator
    %   2- avg_beat_calculator_phase_domain
    %
    % Refs:
    %   1- Sameni, R., Jutten, C., & Shamsollahi, M. B. (2008). Multichannel electrocardiogram 
    %   decomposition using periodic component analysis. 
    %   IEEE Transactions on Biomedical Engineering, 55(8), 1935-1940.
    %   2- Sameni, R., Shamsollahi, M. B., Jutten, C., & Clifford, G. D. (2007). 
    %   A nonlinear Bayesian filtering framework for ECG denoising. IEEE Transactions 
    %   on Biomedical Engineering, 54(12), 2172-2185.
    %   3- B. Vahabzadeh, R. Sameni, The notion of cardiac phase and
    %   its applications in electrophysiological studies, Biomedical Engineering 
    %   (BioMed 2012), Innsbruck, Austria (2012) 110.

    % Authors:
    %   Seyedeh Somayyeh Mousavi 
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA
    
    %%
    if nargin < 4 
            Flag_baseline_alignment = 1;
   end

    try
            % phase calculation
            [phase, ~] = phase_calculator(r_peak_impulse_train);                                      
            % mean ECG extraction
            [feature_vec,~, ~, ~, ~] = avg_beat_calculator_phase_domain(data,phase,num_phase_beans, Flag_baseline_alignment); 
                
            feature_vec = round(feature_vec, 3);
            % Convert Inf to NaN
            feature_vec(isinf(feature_vec)) = NaN;

    catch
            
            feature_vec = nan(1,num_phase_beans);

    end
        
    % Define feature info
    for n = 1:num_phase_beans
                feature_info.names{n} = string(['mean_ecg_beat_sample_', num2str(n)]);
                feature_info.units{n} = "uv";
                feature_info.description{n} = string(['Average ECG beat sample ', num2str(n)]);
    end