from __future__ import division
import os
import numpy as np
import scipy.signal as sig
from .utils import pend_sos, filter_from_npz

__all__ = ['coupling_func', 'witness', 'ideal_estimate']


def coupling_func(true_motion):
    # Bilinear coupling function is just the product of the inputs
    return np.prod(true_motion, axis=0)


def witness(sec, fs, seed=None):
    '''
    Signals meant to represent something like
    beam-spot motion and mirror angulary control signal.
    '''

    if isinstance(seed, np.random.RandomState):
        state = seed
    else:
        state = np.random.RandomState(seed)

    fnyq = fs / 2
    nscale = np.sqrt(fnyq)  # To normalize unit variance to unity ASD

    true_motion = np.zeros((2, fs * sec))
    noise = np.zeros((2, fs * sec))

    # Random noise
    true_motion[0, :] = state.randn(fs * sec) * nscale  # beam spot motion
    true_motion[1, :] = state.randn(fs * sec) * nscale  # ASC control noise

    # Beam spot is driven around by microseism. Realistic spectrum is kind of
    # like f^-2 after microseism until a few Hz, then f^-6 from 3-10 Hz. Noise
    # floor is about 10^-7 lower than ASD at microseism, intercepting at 10 Hz.

    # Make shaping filter for beam spot motion
    # TODO: fix this filter shape to be like ADS error signals
    usos1 = sig.butter(2, [0.1 / fnyq, 0.2 / fnyq],
                       btype='bandpass',
                       output='sos')
    usos2 = sig.butter(4, 3 / fnyq,
                       btype='lowpass',
                       output='sos')
    usos = np.concatenate([usos1, usos2], axis=0)

    # Shape beam spot motion channel
    true_motion[0, :] = sig.sosfilt(usos, true_motion[0, :])

    # Sensing noise of beam spot witness channel
    # Unity ASD to start
    noise[0, :] = state.randn(fs * sec) * nscale
    # Some emperically found scaling to have some SNR over background
    noise[0, :] *= 3e-5

    #  y[0, :] += state.randn()  # Off-center bias
    #  n[0, :] -= state.randn()  # Cancel off-center bias, make it unknown

    # Generate noise with ASC control signal shape
    #true_motion[1, :] = make_ASC_control(sec, fs, seed=state)
    rlp = sig.ellip(4, 3, 40, 17 / fnyq,
                    btype='lowpass',
                    output='sos')

    # We have essentially perfect knowledge of the control force
    # Unity ASD to start
    noise[1, :] = state.randn(fs * sec) * nscale
    # something small like DAC noise
    noise[1, :] *= 1e-13

    # low pass filtering on ASC control signal
    # remove after inserting true CHARD control spectrum
    true_motion[1, :] = sig.sosfilt(rlp, true_motion[1, :])

    # Add sensing noise to witnesses
    witnesses = true_motion + noise

    # Scale empirically chosen to exceed DARM in the tens of Hz
    # not in real units
    # TODO: make this automatically be 10x more than DARM at 10 Hz
    scale = 4e-9
    true_motion *= scale  # scale the angular control noise and spot motion
    witnesses *= scale  # Apply the same scaling to witnesses

    # Turn PUM torque into test mass angle
    true_motion[1, :] = asc_xtal_to_angle(true_motion[1, :], fs)

    return witnesses, true_motion


def ideal_estimate(witnesses, fs):
    '''
    Calculate coupled noise from witness signals, to check regression potential
    '''

    witnesses[1] = asc_xtal_to_angle(witnesses[1], fs)
    estimate = coupling_func(witnesses)

    return estimate


def asc_xtal_to_angle(angular, fs):
    fpend = 8
    fQ = 3
    pend = pend_sos(fpend, fQ, fs)

    # double pendulum
    # ACK! beware of precision noise from too much low pass filtering
    tst_angle = sig.sosfilt(pend, angular)
    tst_angle = sig.sosfilt(pend, tst_angle)

    return tst_angle


def make_ASC_control(sec, fs, seed=None):
    if isinstance(seed, np.random.RandomState):
        state = seed
    else:
        state = np.random.RandomState(seed)

    fnyq = fs / 2
    nscale = np.sqrt(fnyq)  # To normalize unity ASD
    white_noise = state.randn(fs * sec) * nscale

    zpk_file = os.path.join(os.path.dirname(__file__),
                            'ASC_Models/ASC_model.npz')
    sos = filter_from_npz(zpk_file, fs)

    out = sig.sosfilt(sos, white_noise)

    return out
