function extract_ecg_features_multileads(signal_name, data, fs, low_pass_fre, high_pass_fre, ...
                                                                   median_window_coff, mean_window_coff, lead_list, ...
                                                                                                flag_post_processing, output_path)

% function extract_ecg_features_multileads(signal_name, data, fs, low_pass_fre, high_pass_fre, ...
%                                                                  median_window_coff, mean_window_coff, lead_list, ...
%                                                                                               flag_post_processing, output_path)
% Inputs:
%   signal_name (string): Name of the ECG signal
%   data (channels x samples double): raw ECG signal
%   fs (scalar): Sampling frequency of the ECG signal 
%   low_pass_fre (scalar): low-pass cut-off frequency (Hz)
%   high_pass_fre (scalar): high-pass cut-off frequency (Hz)
%   median_window_coff (scalar): coefficient for the moving median window
%   mean_window_coff (scalar): coefficient for the moving mean window
%   lead_list (cell array of strings): Names and order of ECG leads
%   flag_post_processing (binary)(Required for fiducial_det_lsim function)
%   output_path (string): Full path to save the extracted features as a `.csv` file   
%
% Outputs:
%   None. The extracted features are saved directly to the specified .csv file.
%
% Notes:
%   Requires WFDB files in the specified input directory.
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
      [angles_feature_vec, angles_feature_info] = ecg_angle_features(data, fs, ...
                                                                   low_pass_fre, high_pass_fre, ...
                                                                   median_window_coff, mean_window_coff, lead_list, ...
                                                                   flag_post_processing);
      % Define info structs
      feature_infos = { angles_feature_info};
      % Initialize combined metadata
      all_feature_names = {};   
      all_features_units = {};
      all_feature_description = {};

      % Concatenate metadata from each struct
      for i = 1:numel(feature_infos)

               fi = feature_infos{i};
               all_feature_names = [all_feature_names, fi.names];
               all_features_units = [all_features_units, fi.units];
               all_feature_description = [all_feature_description, fi.description];
      end

      % Convert metadata to columns
      feature_names = all_feature_names(:);            
      feature_units = all_features_units(:);           
      feature_description = all_feature_description(:);
      % Convert feature values to a row
      feature_values = num2cell(angles_feature_vec(:)'); 

      % columns 
      feature_data = [
                                feature_description(:)'; 
                                feature_units(:)'; 
                                feature_values
                                ];

      % Row names (vertical)
      row_names = {'Descriptions'; 'Units'; 'Values'};

      % Convert feature names to a row of character vectors
      variable_names = cellstr(feature_names');

      % Table
      Table_csv = cell2table(feature_data, ...
                                 'VariableNames', variable_names');

      % Construct the full file path using output_path and signal_name
      csv_filename = fullfile(output_path, [signal_name '_features_multileads.csv']);

      % Convert the table to a format that can be saved as CSV
      Table_csv = [table(row_names, 'VariableNames', {'Features'}), Table_csv];

      % Write the table to CSV
      writetable(Table_csv, csv_filename);

catch 

      % Error handling
      fprintf("\nERROR while processing\n");
      fprintf("Signal : %s\n", signal_name);

end
end
