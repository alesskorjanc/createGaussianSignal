function plotGaussianSignal(signal, time, mag, freq, typeIn)

% plotGaussianSignal(signal, time, mag, freq, typein)
%
% Functions draws a time-domain (TD) signal generated with function
% gaussianSignal, and a frequency spectrum of the signal defined by the user
% as the gaussianSignal input. It also calculates back the frequency
% spectrum from the output TD signal using Welch's estimator, and draws it
% for comparison.
%
% Created by Ales Skorjanc 18.2.2022
%
%% Inputs:
% 
% signal: row vector
%   TD signal values returned by the function gaussianSignal
%
% time: row vector
%   TD signal time points returned by the function gaussianSignal
%
% mag: row vector
%   amplitude spectral density (asd) or power spectral density (psd) of the
%   signal frequency spectrum defined by the user as gaussianSignal input
%   and returned by the function gaussianSignal in a form interpolated to
%   frequencis in freq
%
% freq: row vector
%   equally spaced frequencies from 0 to half of the sampling frequency
%   (Nyquist frequency) of TD signal, returned by function gaussianSignal
%
% type: string
%   type of the frequency spectrum used as an input into the function
%   gaussianSignal, 'amplitude' or 'power'

%% calculated psd of the TD signal using Welch's estimator

N = length(signal); % length of the TD signal

df = freq(2)-freq(1); % frequency resolution in FD 

fN = max(freq); % Nyquist frequency

fftWindow=2^(nextpow2(N)-0); % fft window length 

averagingWindow=2^(nextpow2(N)-7); % averaging window length

noverlap=0; % averaging window overlap 

% estimate psd of the signal using Welch's estimator, where sampling
% frequency equals 2*fN
[psd, psdf] = pwelch(signal, averagingWindow, noverlap,fftWindow, 2*fN); 

%% recalculate frequency spectrum amplitude from psd and reverse normalization

% if the input into the function gaussianSignal was defined as asd, convert
% psd back to asd by taking the square root of psd and reversing the
% normalization of the signal variance to 1. If the input into
% gaussianSignal was defined as psd, only reverse the normalization of the
% signal variance to 1.

if strcmp(typeIn, 'amplitude')
    magIn = sqrt(psd*sum(mag.^2)*df);  
elseif strcmp(typeIn, 'power')
    magIn = psd*(sum(mag)*df);
end;


%% draw results

figure;

% Plot TD signal
subplot(4, 1, 1);
plot(time, signal);
title('Normalized output signal in time domain')
xlabel('time [s]')
ylabel('amplitude')

% Plot user defined frequency spectrum, used as input for gaussianSignal 
subplot(4, 1, 2);
plot(freq, mag);

if strcmp(typeIn, 'amplitude')
   title('User defined input amplitude spectrum')
   ylabel('asd')
elseif strcmp(typeIn, 'power')
    title('User defined input power spectrum')
    ylabel('psd')
end;

xlabel('frequency [Hz]')


% Plot the amplitude or power spectrum of the TD signal, estimated with
% pwelch, to test the ouput of the gaussianSignal algorithm

subplot(4, 1, 3);
plot(psdf, magIn); 

if strcmp(typeIn, 'amplitude')
   title('Amplitude spectrum of the output TD signal estimated with pwelch')
   ylabel('asd')
elseif strcmp(typeIn, 'power')
    title('Power spectrum of the output TD signal estimated with pwelch')
    ylabel('psd')
end;

xlabel('frequency [Hz]')


% Plot the distribution of values of TD signal to test for their Gaussian
% distribution

subplot(4, 1, 4);
histogram(signal); 
set(gca, 'XLim', [-4, 4]);

title('Distribution of TD signal values')
xlabel('amplitude')
ylabel('N')



