function [data, fs, base_name] = load_and_preprocess_ecg(file_entry, lead_list)

    % function [data, fs, base_name] = load_and_preprocess_ecg(file_entry, lead_list)
    % Load one WFDB ECG record and prepare it for feature extraction 
    %
    % Inputs:
    %   file_entry (struct): one element of a dir()-style file list, with fields .name and .folder
    %   lead_list: expected lead names/order
    %
    % Outputs:
    %   data: preprocessed signal
    %   fs: sampling frequency
    %   base_name: base file name
    %
    % Author: 
    %   Seyedeh Somayyeh Mousavi
    %   Reza Sameni
    %   Emory University, Georgia, USA
    %   bmemousavi@gmail.com
    %
    % AUG, 2026

    %% constant  value
    convert_mv_um = 1000;

    %% Preprocess the data
    % ---- Get file path and base name -----------------------------
    [~, base_name, ~] = fileparts(file_entry.name);
    data_path = file_entry.folder;
    full_file_path = fullfile(data_path, base_name);

    % ---- Load ECG signal -------------------------------------------
    [~, data, fs, siginfo] = rdmat(full_file_path);
    unit = {siginfo.Units};
    gain = {siginfo.Gain};
    challenge_str = {siginfo.Description};

    % Transpose if necessary: signals should be in [channels x samples]
    if size(data, 1) > size(data, 2)
        data = data';
    end
    n_channels = size(data, 1);

    % =========================== Unit Conversion =========================
    for h = 1:n_channels

                unit_lower = lower(unit{h});
                if strcmp(unit_lower, 'mv')
                    convert_ratio = convert_mv_um;
                elseif strcmp(unit_lower, 'uv')
                    convert_ratio = 1;
                else
                    error("Unknown unit in signal %s (channel %d): '%s'", base_name, h, unit{h});
                end
                data(h, :) = convert_ratio * data(h, :);
    end

    % =========================== Gain and Unit Consistency Check ============
    if ~all(strcmp(unit, unit{1})) || ~all(cellfun(@(x) isequal(x, gain{1}), gain))
        error("Units or gains are not consistent in file: %s", base_name); 
    end
    % =========================== Lead Order Check & Reorder ===============
    if ~isequal(lower(challenge_str), lower(cellstr(lead_list)))
            
        [~, new_order] = ismember(lower(lead_list), lower(challenge_str));
        if any(new_order == 0)
            error("Mismatch in channel order or missing leads in: %s", base_name);
        end
        data = data(new_order, :);  
    end

end