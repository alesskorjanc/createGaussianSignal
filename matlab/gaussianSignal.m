function [signalOut, timeOut, magOut, freqOut, typeIn] = gaussianSignal(typeIn, magIn, freqIn, duration, fs)

% [signalOut, timeOut, magOut, freqOut, typeIn] = gaussianSignal(typeIn, magIn, freqIn, duration, fs)
%
% 
% Function returns a time-domain (TD) signal with user-defined spectral
% properties in frequency-domain (FD), using inverse fast Fourier transform
% (ifft).
%
% Created by Ales Skorjanc 18.2.2022
%
%% Inputs:
% 
% typeIn: string
%   user defined type of the input frequency spectrum, 'amplitude' or
%   'power'
%
% magIn: row vector
%   user defined amplitude of the signal frequency spectrum, expressed as
%   amplitude spectral density (asd) per unit bandwidth (1/N), if type
%   'amplitude' is chosen; or signal power spectral density (psd) per
%   frequency bandwidth (Hz) if type 'power' is chosen
%
% freqIn: row vector
%   frequencies at which the signal asd or psd values are defined. Vectors
%   magIn and freqIn must be of the same length and the values in freqIn
%   must increase monotonically.
%
% duration: double
%   duration of the output TD signal in seconds
%
% fs: double
%   sampling frequency of the output TD signal in Hz 
%
%
%% Outputs: 
%
% signalOut: row vector
%   output TD signal with the standard deviation of the amplitude
%   normalized to 1
%
% timeOut: row vector
%   time points at which TD signal was calculated
%
% magOut: row vector
%   user defined input asd or psd, interpolated to frequency values in
%   freqOut
%
% freqOut: row vector
%   equally spaced frequencies from 0 to fs/2 Hz
%
% typeIn: string
%   type of the input frequency spectrum, 'amplitude' or 'power'
% 
% Outputs magOut, freqOut and typeIn can be used with function
% drawGaussianSignal() to visualize the input frequency spectrum and the
% results of the function gaussianSignal.

%% generate the TD signal using ifft

N = round(duration*fs); % length of the output TD signal
fN = fs/2; % Nyquist frequency
dt = 1/fs; % time resolution of the signal in TD

fftL = round(N/2); % length of the signal in FD at positive frequencies 
df = fN/fftL;      % frequency resolution in the FD 

freqOut = 0:df:(fN-df); % positive frequencies in FD from 0 to fN in Hz 

% user defined signal asd or psd interpolated to frequencies in freqOut
magOut =  interp1(freqIn, magIn, freqOut); 

% substitute NaN values, which interp1 returns for points lying outside the
% frequency range defined by the user
magOut(isnan(magOut)) = 0;  

% if user defines asd, convert to power ()^2 and normalize, such that the
% integral of psd over the frequency range 0 to fN Hz equals 1, giving TD
% signal variance 1. If user defines psd, only normalize. 
if strcmp(typeIn, 'amplitude')
    fftP = magOut.^2/(sum(magOut.^2)*df);   
elseif strcmp(typeIn, 'power')
    fftP = magOut/(sum(magOut)*df);
end;

% transform psd to asd, which is the input into the ifft algorithm. Factor
% N comes from Parseval's theorem to ensure that the integral over the
% power spectrum equals the energy of the signal. Factor 2 is needed to
% take into account the doubling of the number of value, which is done below.
fftA = sqrt(fftP*df/2)*N; 

phase = rand(1, fftL)*2*pi(); % generate random phases of FD signal

% calculate complex FD signal with user defined asd or psd and random
% phases
signalFD = fftA.*exp(1i*phase); 

signalFDneg = fliplr(conj(signalFD)); % FD signal values at negative frequencies

signalFD = [signalFD, signalFDneg];    % concatenate positive and negative frequencies

signalFD = [0, signalFD]; % add zero to the beginning of the complex signal

% calculate output TD signal from complex values of the signal in FD using
% ifft
signalOut = ifft(signalFD); 

% generate row vector of times at which the TD signal is defined
timeOut = 0:dt:dt*(length(signalOut)-1); 

% plot the results
plotGaussianSignal(signalOut, timeOut, magOut, freqOut, typeIn)
