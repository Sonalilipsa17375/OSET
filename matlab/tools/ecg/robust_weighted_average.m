function [mn, vr_mn, varargout] = robust_weighted_average(x)
% robust_weighted_average - Robust weighted averaging
%
% Usage:
%   [mn, vr_mn, md, vr_md, best_beat, vr_best_beat, best_beat_index] = robust_weighted_average(x)
%
% Inputs:
%   x: An (N x T) matrix containing N ensembles of a noisy event-related signal of length T
%
% Outputs:
%   mn: The robust weighted average over the N rows of x
%   vr_mn: The variance of the average beat across the N rows of x
%   md: The robust weighted median over the N rows of x (optional)
%   vr_md: The variance of the median beat across the N rows of x (optional)
%   best_beat: The single beat from x closest to the median beat (optional)
%   vr_best_beat: The variance of the best beat across the N rows of x (optional)
%   best_beat_index: index of beat with lowest noise
%
% Reference:
%   J.M. Leski. Robust weighted averaging of biomedical signals. IEEE
%       Trans. Biomed. Eng., 49(8):796-804, 2002.
%
% Revision History:
%   2008: First release
%   2019: Return the average beat variance
%   2021: Return the median beat and its variance
%   2023: Renamed from deprecated version RWAverage
%   2023: Added varargout to speed up when median estimates are not required
%   2026: Added best beat option
%
% Reza Sameni, 2008-2026
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

num_beats = size(x, 1);
if num_beats > 1
    % Average beat (mean)
    mn0 = mean(x, 1);
    noise0 = x - mn0;
    vr = var(noise0, [], 2);
    sm = sum(1 ./ vr);
    weight = 1 ./ (vr * sm);
    mn = weight' * x;
    noise = x - mn;
    vr_mn = var(noise, [], 1);

    if nargout > 2
        % Robust weighted average using the median as the initial estimate
        md0 = median(x, 1);
        noise0 = x - md0;
        vr = var(noise0, [], 2);
        sm = sum(1 ./ vr);
        weight = 1 ./ (vr * sm);
        md = weight' * x;
        varargout{1} = md;
        if nargout > 3
            noise = x - md;
            vr_md = var(noise, [], 1);
            varargout{2} = vr_md;
        end
    end

    if nargout > 4
        % Beat closest to the median beat
        md0 = median(x, 1);
        noise0 = x - md0;
        vr = var(noise0, [], 2);
        [~, selected_beat_index] = min(vr);
        best_beat = x(selected_beat_index, :);
        varargout{3} = best_beat;
        if nargout > 5
            noise = x - best_beat;
            vr_best_beat = var(noise, [], 1);
            varargout{4} = vr_best_beat;
        end
        if nargout > 6
            varargout{5} = selected_beat_index;
        end
    end
else
    mn = x;
    vr_mn = zeros(size(x));
    if nargout > 2
        varargout{1} = x;
    end
    if nargout > 3
        varargout{2} = zeros(size(x));
    end
    if nargout > 4
        varargout{3} = x;
    end
    if nargout > 5
        varargout{4} = zeros(size(x));
    end
    if nargout > 6
        varargout{5} = 1;
    end
end
end