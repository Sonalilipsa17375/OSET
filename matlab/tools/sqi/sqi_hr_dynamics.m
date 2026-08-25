function [rho, mov_rho, sd1, sd2, rr] = sqi_hr_dynamics(peak_indexes, fs, wlen)

rr = diff(peak_indexes)./fs;
wlen_half = round(wlen/2);

% baseline = movmedian(rr, [wlen_half, wlen_half]);
% rr_move_removed = rr - baseline;
% 
% figure
% plot(rr);
% hold on
% plot(baseline);
% grid

rr_move_removed = rr - median(rr);


auto_corr1 = movmean(rr_move_removed(1:end-1).^2, [wlen_half, wlen_half]);
auto_corr2 = movmean(rr_move_removed(2:end).^2, [wlen_half, wlen_half]);
cross_corr = movmean(rr_move_removed(1:end-1) .* rr_move_removed(2:end), [wlen_half, wlen_half]);

mov_rho = cross_corr./sqrt(auto_corr1 .* auto_corr2);
R = corrcoef(rr_move_removed(1:end-1), rr_move_removed(2:end));
rho = R(1, 2);

