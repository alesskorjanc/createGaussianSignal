function plotGaussianSignal(signal,time,mag,freq,typein)

% drawGaussianSignal_Paper_v1(signal,time,mag,freq,typein)
%
% Functions draws the time-domain (TD) signal generated with function
% gaussianSignal,and the frequency spectrum of the signal defined by the
% user as the function input. It also calculates back the frequency
% spectrum of the output TD signal using Welch's estimator, and draws it
% for comparison.
%
% Inputs:
% 
% signal... TD signal returned by function gaussianSignal
%
% time... TD signal time points returned by function gaussianSignal
%
% mag... amplitudes or power spectral density (psd) of the signal frequency
% spectrum defined by the user as gaussianSignal input and returned by the
% function in a form interpolated to frequencis in freq
%
% freq... equally spaced frequencies of the user defined frequency spectrum
% of the signal, returned by function gaussianSignal
%
% type... type of the input frequency spectrum, 'amplitude' or 'power'


%% scale the TD signal to the amplitude used as an input into the function gaussianSignal

N = length(signal);
% N... length 0f the TD signal
df = freq(2)-freq(1);
% df...frequency resolution in FD 
fN = max(freq);
% fN... Nyquist frequency


%% calculated psd of the TD signal using Welch's estimator

fftWindow=2^(nextpow2(N)-4); 
% fftWindow... fft window length 
averagingWindow=2^(nextpow2(N)-7); 
% averagingWindow... averaging window length
noverlap=0; 
% noverlap... window overlap 

[psd,psdf] = pwelch(signal,averagingWindow,noverlap,fftWindow,2*fN); 
% estimate psd of the signal using Welch's estimator

if strcmp(typein,'amplitude')
    magin = sqrt(psd*sum(mag.^2)*df); 
    % if the input into the function gaussianSignal was defined as
    % amplitude spectrum, convert psd back to amplitude by taking the
    % square root of psd and reversing the normalization of the signal
    % variance to 1
elseif strcmp(typein,'power')
    magin = psd*(sum(mag)*df);
    % if the input into the function gaussianSignal was defined as psd,
    % reverse the normalization of the signal variance to 1
end;


%% draw results

figure;

% Plot the generated TD signal
subplot(4,1,1);
plot(time,signal);
title('Normalized output signal in time domain')
xlabel('time [s]')
ylabel('amplitude')

% Plot user defined signal amplitude in FD 
subplot(4,1,2);
plot(freq,mag);

if strcmp(typein,'amplitude')
   title('User defined input amplitude spectrum')
   ylabel('amplitude')
elseif strcmp(typein,'power')
    title('User defined input psd')
    ylabel('psd')
end;

xlabel('frequency [Hz]')


% Plot the amplitude spectrum of the generated signal, estimated with
% pwelch, to test the ouput of the algorithm

subplot(4,1,3);
plot(psdf,magin); 

if strcmp(typein,'amplitude')
   title('Amplitude spectrum of the output TD signal estimated with pwelch')
   ylabel('amplitude')
elseif strcmp(typein,'power')
    title('psd of the output TD signal estimated with pwelch')
    ylabel('psd')
end;

xlabel('frequency [Hz]')


% Plot the distribution of signal's TD amplitudes to test for their
% Gaussian distribution
subplot(4,1,4);
histogram(signal); 
set(gca,'XLim',[-4,4]);
title('Distribution of TD signal amplitudes')
xlabel('amplitude')
ylabel('N')



