#!/usr/bin/env python
from __future__ import division
import sys
import argparse
from mockdata import starting_data, plot_asd, known_models
from scipy.io import savemat
import numpy as np

class helpfulParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = helpfulParser()
parser.add_argument('model', type=str, nargs='?',
                    help="Noise model to use. Defaults to '%(default)s'. "
                         "Valid options are: {}".format(known_models))

parser.add_argument('-t', '--sec', default=16, type=int, nargs='?',
                    help='Integer number of seconds to generate. Defaults '
                         ' to %(default)s')

parser.add_argument('-f', '--fs', default=1024, type=int, nargs='?',
                    help='Sampling frequency of output data in Hz. Defaults '
                         ' to %(default)sHz.')

parser.add_argument('-p', '--doplot', action='store_true', default=False,
                    help='Plot PSDs of background and target.')

parser.add_argument('-k', '--keywords', default=[], nargs=argparse.REMAINDER,
                    help='additional info given to models. Type [model] -k -h '
                         ' for model-specific parser arguments.')

parser.add_argument('-d', '--doshift', action='store_true', default=False,
                    help='Do time shift.')

parser.add_argument('-s', '--shift', default=0.5, type=float, nargs='?',
                    help='Time shift in seconds.')



# Get parameters into global namespace
args   = parser.parse_args()
model  = args.model
sec    = args.sec
fs     = args.fs
doplot = args.doplot
doshift = args.doshift
shift = args.shift
keyword_list = args.keywords

times, background, darm, wit, aux = starting_data(sec=sec, fs=fs,
                                                  model = model,
                                                  parse_str=keyword_list)

if doshift:
    idx = np.mod(np.arange(len(darm),) + int(shift*fs),len(darm))
    darm = darm[idx]

noise_data = {}
noise_data['times']      = times   # just the time array
noise_data['fs']         = fs
noise_data['darm']       = darm     # this is background + nonlin noise
noise_data['wit']        = wit
noise_data['background'] = background
noise_data.update(aux)

# save the dictionary of data into a HDF5 .mat file so that
# its readable in matlab and python

if doshift:
    filename = 'DARM_shift_with_{}'.format(model)
else:
    filename = 'DARM_with_{}'.format(model)
savemat(filename, noise_data,
        appendmat      = True,
        do_compression = True)

if doplot:
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')

    if 'ideal_estimate' in aux:
        est = aux['ideal_estimate']
    else:
        est = None

    plot_asd(background, darm, est, fs=fs, model=model)
