function y = lp_filter_zero_phase_2nd(x, fc)
% lp_filter_zero_phase_2nd - Fourth-order zero-phase Butterworth lowpass filter
%   built from a second-order (biquad) IIR section applied via zero-phase
%   forward-backward filtering.
%
% Syntax: y = lp_filter_zero_phase_2nd(x, fc)
%
% Inputs:
%   x:  Vector or matrix of input data (channels x samples).
%   fc: -3dB cut-off frequency normalized by the sampling frequency
%       (0 < fc < 1). For example, fc = 0.1 corresponds to one tenth
%       of the sampling frequency.
%
% Output:
%   y:  Vector or matrix of filtered data (channels x samples),
%       with zero phase distortion and fourth-order rolloff (-80 dB/decade).
%
% Method:
%   A second-order Butterworth IIR lowpass filter is designed analytically
%   using the bilinear transform (with T=1) and frequency pre-warping at fc.
%   filtfilt() then applies it forward and backward, achieving zero phase
%   and doubling the effective filter order to fourth.
%
%   The analog second-order Butterworth lowpass prototype is:
%
%       H_a(s) = Omega_c^2 / (s^2 + sqrt(2)*Omega_c*s + Omega_c^2)
%
%   The digital cutoff is pre-warped to place the -3dB point exactly at fc:
%
%       Omega_c = 2 * tan(pi * fc)
%
%   Applying the bilinear transform s = 2*(1 - z^{-1})/(1 + z^{-1}) and
%   letting W = Omega_c, a0 = 1 + sqrt(2)*W + W^2, yields:
%
%       b = (W^2 / a0) * [1, 2, 1]
%       a = [1,  2*(W^2 - 1)/a0,  (1 - sqrt(2)*W + W^2)/a0]
%
%   The numerator zeros are at z = -1 (double zero at Nyquist), as expected
%   for a lowpass filter designed via the bilinear transform.
%
%   Revision History:
%       2026: First release, extension of lp_filter_zero_phase()
%
%   Read more:  Mitra, S. (2010). Digital signal processing (4th ed.)
%               New York, NY: McGraw-Hill Professional.
%               Proakis, J. & Manolakis, D. (2006). Digital Signal
%               Processing (4th ed.). Prentice Hall.
%
%   See also: lp_filter_zero_phase, filtfilt, butter
%
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if fc <= 0 || fc >= 1
    error('fc must be in the open interval (0, 1).');
end

% Bilinear transform with frequency pre-warping.
% Maps the digital cutoff fc (normalized, cycles/sample, T=1) to the
% equivalent analog frequency Omega_c so the -3dB point lands exactly
% at fc after discretization.
Wc  = 2 * tan(pi * fc);
Wc2 = Wc^2;

% Second-order Butterworth biquad coefficients derived from the bilinear
% transform of H_a(s) = Wc^2 / (s^2 + sqrt(2)*Wc*s + Wc^2).
% See derivation in the header above.
a0 = 1 + sqrt(2)*Wc + Wc2;

b  = (Wc2 / a0) * [1, 2, 1];       % numerator:   double zero at z = -1 (Nyquist)
a  = [1, ...
      2*(Wc2 - 1) / a0, ...         % a1: real part of pole pair
      (1 - sqrt(2)*Wc + Wc2) / a0]; % a2: conjugate pole product

% Zero-phase filtering via forward-backward application.
% filtfilt squares the magnitude response (4th-order Butterworth rolloff)
% and cancels all phase shift. Input transposed so filtfilt operates along
% the sample dimension (rows); result transposed back to channels x samples.
y = filtfilt(b, a, x')';