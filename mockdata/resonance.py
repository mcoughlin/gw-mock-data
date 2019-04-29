
import numpy as np


def witness(sec=16, fs=2048, w0=2*np.pi*np.sqrt(2)*9.3, Q=100, seed=None):
    '''
    TODO: Add comment

    Parameters
    ----------
    sec: int
        Lengths of time series in seconds.
    fs: int
        Sampling frequency of time series in Hz. Defaults to 2048 Hz
        (2**11)
    w0: int, float
        Resonant frequency. Defaults to 10.
    Q: int, float
        Quality factor. Defaults to 100.
    seed: int, np.random.RandomState
        Random seed

    Returns
    -------
    targets_t:
    witnesses_t:

    '''
    if isinstance(seed, np.random.RandomState):
        state = seed
    else:
        state = np.random.RandomState(seed)

    # Generate random noise
    times = np.linspace(0, sec, fs*sec, endpoint=False)
    dt = times[1] - times[0]  # timestep
    witnesses_t = state.randn(fs * sec)

    # FFT
    freq = np.fft.fftfreq(times.shape[0], d=dt) / (2 * np.pi)
    witnesses_freq = np.fft.fft(witnesses_t, axis=0)

    # Apply transfer function
    targets_freq = witnesses_freq / (w0**2 - freq**2 + w0*freq*1j/Q)

    # inverse FFT
    targets_t = np.fft.ifft(targets_freq)
    targets_t *= 1e-14  # scale

    return targets_t.real, witnesses_t.real
