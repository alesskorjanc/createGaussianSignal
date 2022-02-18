function [signal,time,mag,freq,typein]=gaussianSignal(typein,magin,freqin,duration,fs)

% [signal,time]=gaussianSignal(frequency,amplitude,tmax,fs)
% 
% Function returns a time-domain (TD) signal with user-defined spectral
% properties, using inverse fast Fourier transform (ifft).
%
%% Inputs:
% 
% typein: string
%   user defined type of the input spectrum, 'amplitude' or 'power' 
%
% magin... a row array of user defined amplitudes of the signal frequency
% spectrum, expersed as spectral density per unit bandwidth (1/N), if type
% 'amplitude' is chosen; or signal power spectral density (psd) if type
% 'power' is chosen
%
% freqin... a row array of frequencies at which the signal amplitudes or psd
% are defined. Vectors magin and freqin must be of the same length and
% the values in frequency must increase monotonically.
%
% duration... the duration of the output TD signal in seconds
%
% fs... the sampling frequency of the TD signal in Hz 
%
%
%% Outputs: 
%
% signal... output TD signal with the standard deviation of the amplitude
% normalized to 1.
%
% time... array of TD signal time points
%
% mag... array of user defined input amplitudes or psd of signal
% frequency spectrum, interpolated to frequency values in freq
%
% freq... array of equally spaced frequencies from 0 to fs/2 Hz
%
% type... type of the input frequency spectrum, 'amplitude' or 'power'
% 
% Outputs mag, freq and type can be used with function drawGaussianSignal()
% to visualize the input frequency spectrum and the results.

%% generate the TD signal using ifft

N = round(duration*fs); % length of the output TD signal
fN = fs/2; % Nyquist frequency
dt = 1/fs;  
% dt...time resolution of the signal in TD

fftL = round(N/2);  
% fftL...length of the signal at positive frequencies in FD 
df = fN / fftL; 
% df...frequency resolution in FD 

freq = 0:df:(fN-df);   
% freq...positive frequencies in FD from 0 to fN in Hz 


mag =  interp1(freqin,magin,freq); 
% mag...user defined signal frequency spectrum amplitude or psd
% interpolated to frequencies in freq
mag(isnan(mag)) = 0;  
% substitute NaN values, which interp1 returns for points lying outside the
% frequency range defined by the user


if strcmp(typein,'amplitude')
    % if user defines frequency spectrum amplitude, convert to power and
    % normalize, such that the integral of psd over frequency range equals
    % 1, giving TD signal variance 1
    fftP = mag.^2/(sum(mag.^2)*df); 
elseif strcmp(typein,'power')
    fftP = mag/(sum(mag)*df);
    % if user defines frequency spectrum psd, normalize such that the
    % integral of psd over frequency range equals 1, giving TD signal
    % variance 1
end;

fftA = sqrt(fftP*df/2)*N; 
% transform psd to amplitude density, which is the input into the ifft
% algorithm. Factor N comes from Parseval's theorem to ensure that the
% integral over the power spectrum equals the energy of the signal. Factor 2
% is needed to take into account the doubling of the number of values in
% line 39.

% generate random phases of FD signal:
phase = rand(1,fftL); 

% calculate complex FD signal with user defined amplitudes or psd and
% random phases:
signalFD = fftA.*exp(1i*phase*2*pi()); 

signalFDneg = fliplr(conj(signalFD)); 
% generate FD signal values at negative frequencies
signalFD = [signalFD,signalFDneg]; 
% concatenate positive and negative frequencies
signalFD = [0,signalFD]; 
% add zero to the beginning of the complex signal
signal = ifft(signalFD); 
% calculate output TD signal from complex values of the signal in FD using
% ifft
time = 0:dt:dt*(length(signal)-1); 
% generate vector array of times at which the TD signal is defined


plotGaussianSignal(signal,time,mag,freq,typein)
