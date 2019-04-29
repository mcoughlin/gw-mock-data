from __future__ import division
import numpy as np
import scipy.signal as sig

# Function defaults
sec_d = 16
fs_d = 2048

__all__ = ['coupling_func', 'witness']

# this is the laser wavelength. We need it so that the calibrated witness
# channels make sense; the scattering fringing cares about the ratio of the
# size motion to the wavelength
lambduh = 1064e-9


def coupling_func(y1, y2, scale=3e-18, phi=0):
    return scale*np.sin(2*2*np.pi/lambduh*(y1 + y2) + phi)


def witness(n_relevant, n_irrelevant, sec=sec_d, fs=fs_d, seed=None, rand_phase=False,
            filt_seismic=False):
    '''
    Deterministically generate mock witness channel data for use
    in testing subtraction algorithms.
    '''
    if isinstance(seed, np.random.RandomState):
        state = seed
    else:
        state = np.random.RandomState(seed)

    nyq = fs / 2.0
    times = np.linspace(0, sec, fs*sec, endpoint=False)

    # Let's make this the 'seismic' channel with the low frequency stuff
    y1     = state.randn(sec*fs)
    b1, a1 = sig.butter(2, [0.1 / nyq, 5 / nyq], btype='bandpass')
    y1     = sig.lfilter(b1, a1, y1)
    y1    *= 5e-7    # this puts it into units of meters
    y1_imit = 5e-7*sig.lfilter(b1, a1, state.randn(sec*fs))
    if(filt_seismic):
        b_lp, a_lp = sig.butter(1, 1.1/nyq, btype='lowpass')
        y1_given = sig.lfilter(b_lp, a_lp, y1)
        y1_irr = sig.lfilter(b_lp, a_lp, y1_imit)
    else:
        y1_given = y1
        y1_irr = y1_imit

    # make several witnesses as requested
    rel_1 = (n_relevant+1)//2
    irrel_1 = (n_irrelevant+1)//2
    rel_1_wits = np.tile(y1_given, (rel_1, 1))
    irrel_1_wits = np.tile(y1_irr, (irrel_1, 1))
    wits_1 = np.concatenate((rel_1_wits, irrel_1_wits),0)
    # add different random noise to witnesses
    noise_add1 = 0.0*5e-11*state.randn(rel_1+irrel_1, sec*fs)
    noise_add1 *= state.uniform(low=0.1, high=1.0, size=(rel_1+irrel_1, 1))
    wits_1 += noise_add1

    # this is the acoustic channel so it has only high frequency noise
    f1 = 59.5
    f2 = 119.7
    phase1 = 0.8
    phase2 = 0.32
    y2 =  1   * np.sin(2*np.pi * (f1 * times + phase1))  # a couple of sine waves
    y2 += 0.3 * np.sin(2*np.pi * (f2 * times + phase2))
    y2 *= 1e-7       #  this puts it into units of meters

    rel_2 = n_relevant//2
    irrel_2 = n_irrelevant//2
    if(rand_phase):
        phase1_wit = 2*np.pi*state.uniform(size=(rel_2, 1))
        phase2_wit = 2*np.pi*state.uniform(size=(rel_2, 1))
    else:
        phase1_wit = np.zeros((rel_2, 1))
        phase2_wit = np.zeros((rel_2, 1))

    rel_2_wits = 1e-7*np.sin(2*np.pi*f1*times + phase1 + phase1_wit) + \
                 3e-8*np.sin(2*np.pi*f2*times + phase2 + phase2_wit)

    f1_irr = 71.2
    f2_irr = 143.0
    phase1_irr = 2*np.pi*state.uniform(size=(irrel_2, 1))
    phase2_irr = 2*np.pi*state.uniform(size=(irrel_2, 1))
    irrel_2_wits = 1e-7*np.sin(2*np.pi*f1_irr*times + phase1_irr) + \
                   3e-8*np.sin(2*np.pi*f2_irr*times + phase2_irr)
    wits_2 = np.concatenate((rel_2_wits, irrel_2_wits), 0)
    wits = np.concatenate((wits_1, wits_2), 0)

    return y1, y2, wits
