from __future__ import division
import numpy as np

from .mock_bg import bucket_noise, get_lines
from . import scatter
from . import bilinear
from . import resonance
import argparse

__all__ = ['starting_data', 'known_models']

known_models = [
    'scatter',
    'bilinear',
    'resonance',
]

# Function defaults
sec_d = 16     # default duration in seconds
fs_d  = 2048   # sample rate for fast channels (DARM, ASC, slow)
asd_d = 2e-20  # noise level (m/rHz) of bucket


def starting_data(sec=sec_d, fs=fs_d, model='scatter', seed=None,
                  parse_str=''):
    '''
    Integrated function for nonlinear noise subtraction investigations.

    Generates DARM background, and adds nonlinear noise of a given
    model.

    Parameters
    ----------
    sec : integer
        Length of time series in seconds. Defaults to 16
    fs : integer, optional
        Sampling frequency of time series in Hz. Defaults to 2048 Hz
        (2**11)
    model : string, optional
        String that dictates the nonlinear noise model being generated. See
        `Models` section below for currently implemented models. Defaults to
        'scatter'.
    seed : int or np.random.RandomState instance, optional
        If an integer, used as the seed for a new
        `np.random.RandomState` instance that generates the timeseries.
        If an instance of `RandomState`, it is used directly. If `seed`
        is ``None``, the `RandomState` will try to read data from
        ``/dev/urandom`` (or the Windows analogue) if available or seed
        from the clock otherwise. Defaults to ``None``.

    Returns
    -------
    times : ndarray
        Array of times
    background : ndarray
        Array of mock DARM background with no additional noise added, to
        be used to evaluate the regression efficacy.
    target : ndarray
        Background with nonlinear noise added to be used as the
        subtraction target
    wit : ndarray, shape (N, tar.size)
        N witness signals used as input to nonlinear regression methods.
    aux : Dict
        Python dictionary containing auxiliary data that may vary depending on
        the noise model being used. For instance, "perfect" witnesses,
        intermediate signals, or best case regression results.

    Models
    ------

    'scatter' :
        One witness moving slowly, and one acousticly active witness coupling
        together via the form sin(4 pi/lambda*(y1 + y2) + phi)
    '''
    ############################
    # To add a new noise model #
    ############################
    #  - Add the string ID to the `known_models` list at the top of the file
    #  - Add an `elif` to the if block and call your model-specific functions
    #    in there. I.e.:
    #      - Add any keywords or aruments your model needs to the parser
    #      - Generate witness data streams
    #      - Define and apply any nonlinear functions & filters
    #      - Create the list of witnesses to be used in the regression in a
    #        list called `witnesses`
    #      - Create the coupled noise time series in `coupled_noise`
    #      - Add any additional data you want out of the function to `aux`

    if model not in known_models:
        raise ValueError('Unknown noise model type: {}'.format(model))

    if isinstance(seed, np.random.RandomState):
        state = seed
    else:
        state = np.random.RandomState(seed)

    times = np.linspace(0, sec, fs * sec, endpoint=False)
    background = bucket_noise(sec, fs, seed=state)
    background += get_lines(sec, fs, seed=state)

    aux = {}  # Define in case we don't want to use it.
    parser = argparse.ArgumentParser(prog='-k')

    if model == 'scatter':
        parser.add_argument(
            '-p',
            '--random_phase',
            action='store_true',
            help='Invoke to randomize phase of audio band sine waves.')
        parser.add_argument(
            '-f',
            '--filt_seismic',
            action='store_true',
            help='Invoke to apply low-pass filter to seismic witness.')
        parser.add_argument(
            '-r',
            '--relevant',
            type=int,
            default=2,
            help='Number of relevant witnesses to include.')
        parser.add_argument(
            '-i',
            '--irrelevant',
            type=int,
            default=0,
            help='Number of irrelevant witnesses to include.')

        args = parser.parse_args(parse_str)
        rand_phase = args.random_phase
        filt_seismic = args.filt_seismic
        rel = args.relevant
        irrel = args.irrelevant

        y1, y2, wits = scatter.witness(rel, irrel,
                                       sec=sec,
                                       fs=fs,
                                       seed=state,
                                       rand_phase=rand_phase,
                                       filt_seismic=filt_seismic)

        coupled_noise = scatter.coupling_func(y1, y2)

        witnesses = wits
        aux['y1'] = y1
        aux['y2'] = y2

    elif model == 'bilinear':
        parser.add_argument(
            '-p',
            '--pairs',
            type=int,
            default=1,
            help='Number of beam spot + angular motion channel pairs to '
                 'return, all of which contribute noise via the bilinear '
                 'coupling.')

        args = parser.parse_args(parse_str)
        pairs = args.pairs

        if pairs < 1 or pairs % 1 != 0 :
            raise ValueError('Pairs must be a positive integer')

        coupled_noise = 0
        ideal_estimate = 0
        beam_motions = []
        angular_controls = []
        true_beam = []
        true_angular = []

        for _ in range(pairs):
            witness, true_motion = bilinear.witness(sec, fs, seed=state)

            coupled_noise += bilinear.coupling_func(true_motion)
            ideal_estimate += bilinear.ideal_estimate(witness, fs)

            beam_motions.append(witness[0])
            angular_controls.append(witness[1])

            true_beam.append(true_motion[0])
            true_angular.append(true_motion[1])

        witnesses = np.stack(beam_motions + angular_controls)
        true_motions = np.stack(true_beam + true_angular)

        aux = {}
        aux['true_motions'] = true_motions
        aux['ideal_estimate'] = ideal_estimate
        aux['Npairs'] = pairs

    elif model == 'resonance':
        default_freq = 9.3*np.sqrt(2)
        parser.add_argument('-q','--quality',
                            type=float,
                            dest='quality',
                            default=100,
                            help='quality factor')
        parser.add_argument('-f', '-frequency',
                            type=float,
                            dest='frequency',
                            default=default_freq,
                            help='resonant frequency (Hz)')
        args = parser.parse_args(parse_str)

        print('Quality: {}'.format(args.quality))
        print('Resonant frequency: {}'.format(args.frequency))

        w0 = 2*np.pi*args.frequency
        coupled_noise, witnesses = resonance.witness(
            sec=sec, fs=fs, w0=w0, Q=args.quality)

    target = coupled_noise + background

    return times, background, target, witnesses, aux
