from __future__ import division
import numpy as np
import scipy.signal as sig

# These are the things that get imported when running `from foo import *`
__all__ = ['bucket_noise', 'get_lines']


def bucket_noise(
        sec,
        fs,
        norm_freq=100,
        norm_amp=2.2e-20,
        zeros=None,
        poles=None,
        seed=None, ):
    """
    Generate noise with approximate spectral properties of aLIGO noise
    curve.

    Parameters
    ----------
    sec : integer
        Length of time series in seconds.
    fs : integer
        Sampling frequency of time series in Hz. Defaults to 2048 Hz
    norm_freq : float
        Frequency where ASD normalization is specified in Hz. Defaults
        to 100 Hz.
    norm_amp : float
        Desired amplitude spectral density at `norm_freq` in 1/sqrt(Hz).
        Defaults to 4e-20 1/sqrt(Hz).
    zeros : array_like
        List or array of zeros that describe the background ASD shape,
        in Hz. If `None`, uses 8 zeros at 20 Hz and one zero at 350 Hz.
        Defaults to `None`.
    poles : array_like
        List or array of poles that describe the background ASD shape,
        in Hz. If `None`, uses 8 poles at 9 Hz. Defaults to `None`.
    seed : int or np.random.RandomState instance, optional
        If an integer, used as the seed for a new
        `np.random.RandomState` instance that generates the timeseries.
        If an instance of `RandomState`, it is used directly. If `seed`
        is ``None``, the `RandomState` will try to read data from
        ``/dev/urandom`` (or the Windows analogue) if available or seed
        from the clock otherwise. Defaults to ``None``.

    Returns
    -------
    data : ndarray, shape (sec*fs,)
        Time series signal imitating LIGO noise.

    """

    if isinstance(seed, np.random.RandomState):
        state = seed
    else:
        state = np.random.RandomState(seed)

    if poles is None:
        poles = [5] * 2
    if zeros is None:
        zeros = [24] * 2 + [350]

    N = int(sec * fs)
    Nfft = _nextpow2(N)

    nscale = fs * np.sqrt(Nfft / fs / 2)  # To normalize to unity ASD

    angles = state.uniform(low=-np.pi, high=np.pi, size=Nfft // 2)
    rfft_data = np.zeros(Nfft // 2 + 1, dtype=np.complex128)
    rfft_data[1:] = nscale * np.exp(1j * angles)

    freqs = np.fft.rfftfreq(Nfft, d=1 / fs)
    shape = _asd_shape(freqs, zeros, poles, norm_freq, norm_amp)

    data = np.fft.irfft(rfft_data * shape)

    return data[:N]


def _asd_shape(freqs, zeros, poles, norm_freq, norm_amp):

    poles = -1j * np.asarray(poles)[:, np.newaxis]
    zeros = -1j * np.asarray(zeros)[:, np.newaxis]

    norm = norm_amp * np.prod(norm_freq - poles) / np.prod(norm_freq - zeros)
    resp = np.prod(freqs - zeros, axis=0) / np.prod(freqs - poles, axis=0)
    return norm * resp


def _nextpow2(n):
    p = int(np.ceil(np.log2(n)))
    return 2**p


def get_lines(sec, fs, peak_amp=1e-19, seed=None):
    """
    Generate line noise at 60Hz and harmonics

    Parameters
    ----------
    sec : integer
        Length of time series in seconds.
    fs : integer
        Sampling frequency of time series in Hz. Defaults to 2048 Hz
    peak_amp : float
        Amplitude of strongest harmonic. Defaults to 1e-19.
    seed : int or np.random.RandomState instance, optional
        If an integer, used as the seed for a new
        `np.random.RandomState` instance that generates the timeseries.
        If an instance of `RandomState`, it is used directly. If `seed`
        is ``None``, the `RandomState` will try to read data from
        ``/dev/urandom`` (or the Windows analogue) if available or seed
        from the clock otherwise. Defaults to ``None``.

    Returns
    -------
    data : ndarray, shape (sec*fs,)
        Time series signal
    """

    if isinstance(seed, np.random.RandomState):
        state = seed
    else:
        state = np.random.RandomState(seed)

    line_freqs = [60, 120, 180]
    relative_amplitudes = [1, .2, .5]
    line_phases = state.uniform(-np.pi, np.pi, size=len(line_freqs))

    maxA = 0
    tt = np.arange(sec * fs) / fs
    out = np.zeros_like(tt)

    for ii, freq in enumerate(line_freqs):
        if freq < 0.8 * (fs / 2):  # Sub-nyquist, ok to add line

            A = relative_amplitudes[ii]
            w = 2 * np.pi * freq
            phi = line_phases[ii]
            out += A * np.sin(w * tt + phi)
            maxA = max(maxA, A)

    if maxA > 0 :
        out *= peak_amp / maxA

    return out


def iir_bp(fstops, fs, order=4):
    """
    the opposite of a notch filter: return coefficients to add a spectral line
    """
    nyq = 0.5 * fs

    # Zeros zd, poles pd, and gain kd for the digital filter
    zd = np.array([])
    pd = np.array([])

    # Lines
    for fstopData in fstops:
        fstop = fstopData[0]
        df = fstopData[1]
        df2 = fstopData[2]
        low = (fstop - df) / nyq
        high = (fstop + df) / nyq
        low2 = (fstop - df2) / nyq
        high2 = (fstop + df2) / nyq
        z, p, k = sig.iirdesign([low2, high2], [low, high],
                                gpass=1,
                                gstop=6,
                                ftype='ellip',
                                output='zpk')
        zd = np.append(zd, z)
        pd = np.append(pd, p)

    # Return the numerator and denominator of the digital filter
    b, a = sig.zpk2tf(zd, pd, k)
    return b, a
