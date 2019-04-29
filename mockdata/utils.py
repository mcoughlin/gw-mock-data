from __future__ import division
#  import os
import numpy as np
import scipy.signal as sig


def filter_from_npz(filename, fs):
    with np.load(filename) as data:
        z = data['z']
        p = data['p']
        k = data['k']

    #sos conversion
    zd, pd, kd = zpk_bilinear(z, p, k, fs)
    sys_disc = sig.ZerosPolesGain(zd, pd, kd, dt=1 / fs)
    sos = sig.zpk2sos(sys_disc.zeros, sys_disc.poles, sys_disc.gain)
    return sos


def zpk_bilinear(z, p, k, fs, f_warp=None):
    '''
    Perform discretization of a continuous ZPK system.

    Uses bilinear transform, optionally prewarping at `f_warp`.
    '''
    # Manually do bilinear transform
    # If you are troubled by this, check out the MATLAB source for the
    # `bilinear` function at matlab/toolbox/signal/signal/bilinear.m

    if f_warp is not None:
        w_warp = 2 * np.pi * f_warp
        f2 = w_warp / np.tan(w_warp / (fs * 2))
    else:
        f2 = 2 * fs

    z = z[np.isfinite(z)]  # Remove zeros at infinity

    zd = (1 + z / f2) / (1 - z / f2)
    pd = (1 + p / f2) / (1 - p / f2)
    kd = np.real(k * np.prod(f2 - z) / np.prod(f2 - p))

    # Add zeros at z=-1 to get same number of poles and zeros
    zd = np.concatenate((zd, -1*np.ones(pd.size-zd.size)))

    return zd, pd, kd


def pend_sos(f0, Q, fs, dc_gain = 1):
    '''
    Make a digital filter for a pendulum TF in SOS form.
    '''

    if f0 >= fs / 2:
        raise ValueError('f0 too high for nyquist')
    if Q < 0:
        raise ValueError('Q must be >=0')

    w_pole = -2 * np.pi * f0
    angle = np.arccos(1 / 2 / Q)

    zeros = np.array([])
    poles = w_pole * np.array([np.exp(1j * angle), np.exp(-1j * angle)])
    gain = dc_gain * np.prod(poles)

    #sos conversion
    zd, pd, kd = zpk_bilinear(zeros, poles, gain, fs)
    sys_disc = sig.ZerosPolesGain(zd, pd, kd, dt=1 / fs)
    sos = sig.zpk2sos(sys_disc.zeros, sys_disc.poles, sys_disc.gain)

    return sos

