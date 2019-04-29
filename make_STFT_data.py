import noisesub as n
from mockdata import starting_data
import numpy as np
from regplots import plot_results

times, data, darm, wit, _ = starting_data(model='scatter', sec=2048)
wit1 = wit[0]
wit2 = wit[1]

print "starting data length", darm.shape, times.shape
fs = int(1/times[1])

plot_results(darm, wit1, wit2, darm-data, data, plot_all=False)
freqs, ts, STFT_darm = n.stft(darm, fs=fs, nperseg=4096,
    return_onesided=True)

freqs, ts, STFT_wit1 = n.stft(wit1, fs=fs, nperseg=4096,
    return_onesided=True)

freqs, ts, STFT_wit2 = n.stft(wit2, fs=fs, nperseg=4096,
    return_onesided=True)

freqs, ts, STFT_data = n.stft(data, fs=fs, nperseg=4096,
    return_onesided=True)


np.save("STFT_wit1", STFT_wit1)
np.save("STFT_wit2", STFT_wit2)
np.save("STFT_darm", STFT_darm)
np.save("STFT_data", STFT_data)
