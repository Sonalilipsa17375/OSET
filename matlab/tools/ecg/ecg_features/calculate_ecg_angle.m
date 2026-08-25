function angle = calculate_ecg_angle(A, B)

    % function angle = calculate_ecg_angle(A, B)   
    % Formula: atan2(A, B) (degree)
    % Computes the angle (in degrees) between two ECG signal components
    %
    % Inputs:
    %   A: numerator component
    %   B: denominator component
    %
    % Outputs:
    %   angle in degrees
    %
    % Authors:
    %   Seyedeh Somayyeh Mousavi 
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA

    if isfinite(A) && isfinite(B)
        angle = rad2deg(atan2(A, B));
    else
        angle = NaN;
    end

end