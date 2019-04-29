#!/usr/bin/env python

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from mockdata.mock_bg import bucket_noise, add_lines

fs = 4096*4
data = bucket_noise(16, fs)
data = add_lines(data, fs)
NFFT = fs/2
# make ASD plot to verify noise spectrum
Pxx_data, freqs = mlab.psd(data, Fs = fs, NFFT = NFFT)
asd_data = np.sqrt(Pxx_data)

plt.figure()
plt.loglog(freqs,  asd_data,'b',label='mock noise background')
plt.xlim([10, fs/2])
plt.grid('on')
plt.ylabel('ASD (strain/rtHz)')
plt.legend(loc='upper center')
plt.xlabel('Freq (Hz)')
#plt.savefig('noisecurve.png')
plt.show()
