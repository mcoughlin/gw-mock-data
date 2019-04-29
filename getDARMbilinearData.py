#!/usr/bin/env python
# this function gets some data (from the 40m) and saves it as
# a .mat file for the matlabs
# Ex. python -O getData.py

from __future__ import division

import os, sys, time
import numpy as np

import scipy.io as sio
import scipy.signal as sig
from astropy.time import Time

nds_osx = '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/'
nds_sandbox = '/usr/lib/python2.7/dist-packages/'
if os.path.exists(nds_osx):
    sys.path.append(nds_osx)
elif os.path.exists(nds_sandbox):
    sys.path.append(nds_sandbox)

import nds2

if isinstance(sys.argv[1], str):
    ifo = sys.argv[1]
else:
    sys.exit("input argument must be a string")

doPlots = False
if doPlots:
    import gwpy.time, gwpy.timeseries
    import gwpy.frequencyseries, gwpy.spectrogram
    import gwpy.plotter
 
    plotDirectory = 'plots'
    if not os.path.isdir(plotDirectory):
        os.makedirs(plotDirectory)

# channel names
chan_head = ifo + ':'
#fname = 'ChanList_S2L.txt'
fname = 'ChanList_darm.txt'

with open(fname, 'r') as f:
    chanlines = f.read().split()
channels = [chan_head + line for line in chanlines]


# good times / bad times
# times      = '2017-01-12 08:00:00'   # this is early January when ASC coupling was high

# 3 hour stretch with moderate high uSeism
if ifo == 'L1':
    ndsServer  = 'nds.ligo-la.caltech.edu'
    times      = '2017-01-04 11:40:00'
elif ifo == 'H1':
    ndsServer  = 'nds.ligo-wa.caltech.edu'
    times      = '2017-01-04 11:40:00'
    #times      = '2017-01-04 10:05:00'

else:
    sys.exit("unknown IFO specified")

# Setup connection to the NDS
portNumber = 31200
conn       = nds2.connection(ndsServer, portNumber)

# Setup start and stop times
# good double coinc time = '2017-03-11 14:00:00'
#times   = '2017-03-11 14:00:00'
t       = Time(times, format='iso', scale='utc')
t_start = int(t.gps)
dur     = 2048

# Data will be downsampled to `fsup` Hz
fsup = 256
if __debug__:
    print("Output sample rate: {} Hz".format(fsup))


if __debug__:
    print("List of Channels:...")
    print("\n".join(channels))

print("Getting data from " + ndsServer + "...")
tic = time.time()
data = conn.fetch(t_start, t_start + dur, channels)

# get the data and stack it into a single matrix where the data are the columns
vdata = []
for k in range(len(channels)):
    fsdown = data[k].channel.sample_rate
    down_factor = int(fsdown // fsup)

    fir_aa = sig.firwin(20 * down_factor + 1, 0.8 / down_factor,
        window='blackmanharris')

    # Prevent ringing from DC offset
    DC = np.mean(data[k].data)

    # Using fir_aa[1:-1] cuts off a leading and trailing zero
    downdata = sig.decimate(data[k].data, down_factor,
                            ftype = sig.dlti(fir_aa[1:-1], 1.0),
                            zero_phase = True)
    vdata.append(downdata+DC)

    if doPlots:
        pngFile = os.path.join(plotDirectory,"%s.png"%(channels[k].replace(":","_")))
        dataTime = gwpy.timeseries.TimeSeries(downdata, sample_rate = fsdown, epoch = t_start, dtype=float)
        NFFT = 8.0
        overlap = 4.0
        dataASD = dataTime.asd(fftlength=NFFT,overlap=overlap,method='welch')

        plot = dataASD.plot()
        ax = plot.gca()
        plot.xlim = [10.0,256]
        #plot.ylim = [10**-10, 10**-4]
        plot.xlabel = "Frequency [Hz]"
        plot.axes[0].set_xscale("log")
        plot.axes[0].set_yscale("log")
        plot.save(pngFile,dpi=200)
        plot.close()

# save to a hdf5 format that matlab can read (why is compression off by default?)
funame = 'Data/' + ifo + '_data_array.mat'
#funame = 'Data/' + ifo + '_darm.mat'

sio.savemat(funame,
            mdict={'data': vdata, 'fsample': fsup, 'chans': channels},
            do_compression=True)
print("Data saved as " + funame)

if __debug__:
    print("Channel name is " +          data[0].channel.name)
    print("Sample rate is " +       str(data[0].channel.sample_rate) + " Hz")
    print("Number of samples is " + str(data[0].length))
    print("GPS Start time is " +    str(data[0].gps_seconds))
    print("Data retrieval time = " + str(round(time.time() - tic,3)) + " s")

# uncomment this stuff to get info on what fields are in the data
#dir(data[0])
#dir(data[0].channel)
