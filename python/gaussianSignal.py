import numpy as np
from scipy.signal import welch
import matplotlib.pyplot as plt


def gaussianSignal(spec_type, spec_mag, spec_freq, duration, fs):
    """Return a Gaussian noise with user-defined amplitude spectrum.

    Created by Jan Benda 12.10.2022 based on the matlab functions by Ales Skorjanc.

    Parameters
    ----------
    spec_type: str
        Type of the input frequency spectrum `spec_mag`: 'amplitude' or 'power'.
    spec_mag: array of floats
        Amplitude of the frequency spectrum, expressed as amplitude
        spectral density (asd) per unit bandwidth (1/N), if
        `spec_type` is 'amplitude', or power spectral density (psd)
        per frequency bandwidth (Hz) if `spec_type` is 'power'.
    spec_freq: array of floats
        Frequencies at which the asd or psd values in `spec_mag` are
        defined. Arrays `spec_mag` and `spec_freq` must be of the same
        length and the values in `spec_freq` must increase
        monotonically, but do not need to be equidistant.
    duration: float
        Duration of the generated Gaussian noise in seconds.
    fs: floate
        Sampling frequency of the generated Gaussian noise in Hz.

    Returns
    -------
    noise_sig: array of floats
        The generated Gaussian noise signal with standard deviation
        normalized to 1
    times: array of floats
        Time points at which the Gaussian noise `noise_sig` was calculated.
    spec_mag: array of floats
        The user supplied amplitude or powe spectrum, interpolated to
        the frequency values in `spec_freq`.
    freqs: array
        Equally spaced frequencies for `spec_mag` from 0 to fs/2 Hz
    spec_type: str
        Type of the frequency spectra: 'amplitude' or 'power'
    
    Outputs spec_mag, freqs and spec_type can be used with function
    drawGaussianSignal() to visualize the input frequency spectrum and the
    results of the function gaussianSignal.
    """

    # the time domain (TD) signal using inverse fast Fourier Transform (FFT):
    
    N = np.round(duration*fs) # length of the output TD signal
    fN = fs/2                 # Nyquist frequency
    dt = 1/fs                 # time resolution of the signal in TD

    fftL = round(N/2)         # length of the signal in frequency domain (FD) at positive frequencies 
    df = fN/fftL              # frequency resolution in the FD 

    freqs = np.arange(0, fN, df) # positive frequencies in FD from 0 to fN in Hz

    # user defined spectrum interpolated to frequencies in `freqs`:
    spec_mag =  np.interp(freqs, spec_freq, spec_mag, left=0, right=0)
    spec_mag[0] = 0           # no offset

    # if user defines asd, convert to power ()^2 and normalize, such that the
    # integral of psd over the frequency range 0 to fN Hz equals 1, giving TD
    # signal variance 1. If user defines psd, only normalize. 
    if spec_type == 'amplitude':
        fftP = spec_mag**2/np.sum(spec_mag**2)/df
    elif spec_type == 'power':
        fftP = spec_mag/np.sum(spec_mag)/df
    else:
        raise ValueError(f'Invalid value "{spec_type}" for spec_type. Must be either "amlpitude" or "power".')

    # transform psd to asd, which is the input into the ifft algorithm. Factor
    # N comes from Parseval's theorem to ensure that the integral over the
    # power spectrum equals the energy of the signal. Factor 2 is needed to
    # take into account the doubling of the number of values further below.
    fftA = np.sqrt(fftP*df/2)*N 

    phases = 2*np.pi*np.random.random(fftL)  # generate random phases of FD signal

    # calculate complex FD signal based on amplitude spectrum and random phases:
    signalFD = fftA*np.exp(1j*phases) 

    # transform FD signal back into TD:
    noise_sig = np.fft.irfft(signalFD) 

    # generate times at which the TD signal is defined:
    times = np.arange(len(noise_sig))*dt
    
    return noise_sig, times, spec_mag, freqs, spec_type


def plotGaussianSignal(signal, time, mag, freq, typeIn):
    """Functions draws a time-domain (TD) signal generated with function
    gaussianSignal, and a frequency spectrum of the signal defined by the
    user as the gaussianSignal input. It also calculates back the
    frequency spectrum from the output TD signal using Welch's estimator,
    and draws it for comparison.

    Created by Ales Skorjanc 18.2.2022

    Paramters
    ---------
    signal: row vector
    TD signal values returned by the function gaussianSignal
    
    time: row vector
    TD signal time points returned by the function gaussianSignal
    
    mag: row vector
    amplitude spectral density (asd) or power spectral density (psd) of the
    signal frequency spectrum defined by the user as gaussianSignal input
    and returned by the function gaussianSignal in a form interpolated to
    frequencis in freq
    
    freq: row vector
    equally spaced frequencies from 0 to half of the sampling frequency
    (Nyquist frequency) of TD signal, returned by function gaussianSignal
    
    type: string
    type of the frequency spectrum used as an input into the function
    gaussianSignal, 'amplitude' or 'power'
    """
    # compute psd of the TD signal using Welch's estimator:
    N = len(signal)
    df = freq[1] - freq[0]               # frequency resolution in FD 
    fN = freq[-1]                        # Nyquist frequency
    #fftWindow = 2**nextpow2(N)           # fft window length
    fftWindow = 1024
    noverlap = 0                         # averaging window overlap 
    # estimate psd of the signal using Welch's estimator, where sampling
    # frequency equals 2*fN
    psdf, psd = welch(signal, noverlap=noverlap, nperseg=fftWindow, fs=2*fN) 

    ## recalculate frequency spectrum amplitude from psd and reverse normalization

    # if the input into the function gaussianSignal was defined as asd, convert
    # psd back to asd by taking the square root of psd and reversing the
    # normalization of the signal variance to 1. If the input into
    # gaussianSignal was defined as psd, only reverse the normalization of the
    # signal variance to 1.
    if spec_type == 'amplitude':
        magIn = np.sqrt(psd*np.sum(mag**2)*df)  
    elif spec_type == 'power':
        magIn = psd*np.sum(mag)*df
    else:
        raise ValueError(f'Invalid value "{spec_type}" for spec_type. Must be either "amlpitude" or "power".')

    ## draw results
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    fig, axs = plt.subplots(4, 1, constrained_layout=True)

    # Plot TD signal
    axs[0].plot(time, signal)
    axs[0].set_title('Normalized output signal in time domain')
    axs[0].set_xlabel('time [s]')
    axs[0].set_ylabel('amplitude')

    # Plot user defined frequency spectrum, used as input for gaussianSignal 
    axs[1].plot(freq, mag)
    if spec_type == 'amplitude':
        axs[1].set_title('User defined input amplitude spectrum')
        axs[1].set_ylabel('asd')
    elif spec_type == 'power':
        axs[1].set_title('User defined input power spectrum')
        axs[1].set_ylabel('psd')
    axs[1].set_xlabel('frequency [Hz]')


    # Plot the amplitude or power spectrum of the TD signal, estimated with
    # pwelch, to test the ouput of the gaussianSignal algorithm
    axs[2].plot(psdf, magIn) 
    if spec_type == 'amplitude':
        axs[2].set_title('Amplitude spectrum of the output TD signal estimated with pwelch')
        axs[2].set_ylabel('asd')
    elif spec_type == 'power':
        axs[2].set_title('Power spectrum of the output TD signal estimated with pwelch')
        axs[2].set_ylabel('psd')
    axs[2].set_xlabel('frequency [Hz]')

    # Plot the distribution of values of TD signal to test for their Gaussian
    # distribution
    axs[3].hist(signal, 200)
    axs[3].set_xlim(-4, 4)
    axs[3].set_title('Distribution of TD signal values')
    axs[3].set_xlabel('amplitude')
    axs[3].set_ylabel('N')

    plt.show()


if __name__ == '__main__':
    duration = 10.0
    fs = 10000.0
    spec_freq = np.arange(0, 3000, 0.1)
    spec_mag = 1-np.cos(2*np.pi*spec_freq/1000)
    noise_sig, times, spec_mag, freqs, spec_type = gaussianSignal('amplitude', spec_mag, spec_freq, duration, fs)
    plotGaussianSignal(noise_sig, times, spec_mag, freqs, spec_type)

