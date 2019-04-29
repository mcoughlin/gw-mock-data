import sys
import numpy as np
import scipy.signal as sig
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import noisesub

from mockdata import starting_data, my_func, coupling_filter,
    apply_filt, FS, SEC, NOISE_LEVEL


def coherence_sub(tar, wit, td_seg=True, signif=True):
    '''
    calculate frequency domain coherence between guessed function and target
    '''
    N = tar.size

    bins = 128 # bins per frequency band
    bands = N//(2*bins) # number of frequency bands

    if(td_seg):
        ff, coherence, coefs = noisesub.mcoherence(tar, wit, FS, nperseg=2*bands)
            # noverlap=0, window="boxcar")
        coefs = coefs[0]
    else:
        csd = inner_product(tar, wit, bands)
        tpsd = inner_product(tar, tar, bands)
        wpsd = inner_product(wit, wit, bands)
        coefs = csd/wpsd
        coherence = np.real((csd*np.conj(csd)/(tpsd*wpsd)))

    # if signif=True, don't perform subtraction at freqs with low coherence
    if(signif):
        set_0 = np.less(abs(coherence), .05)
        coefs[set_0] = 0.0
        low = 10*SEC//bins
        coefs[0:low] = np.zeros(low)

    if(td_seg):
        # maybe there's some way to apply windowing to the coefs to match the csd?
        coefs = np.repeat(coefs, bins)
        coefs = coefs[bins//2:-1*bins//2]
        #coefs = np.append(coefs, 0.0)
    else:
        coefs = np.repeat(coefs, bins)
        #coefs = np.append(coefs, 0.0)

    sub = np.fft.rfft(wit)[:-1]
    # discard nyquist bin so that the bands fit evenly
    sub = sub*coefs
    tarfft = np.fft.rfft(tar)[:-1]

    diff = tarfft - sub
    y_guess = np.fft.irfft(diff, len(tar))

    return diff, coherence, y_guess

def coherence_plots(tar, wit, signif, file_end=""):
    '''
    Generate and save plots showing:
    (1) coherence between estimated noise and original strain
        data as a function of frequency. Coherence below 0.05 is not considered
        significant.
    (2) estimated magnitude of transfer function applied to noise before
        it was added to the background, compared with transfer function of
        actual filter applied.
    '''
    N = tar.size

    bins = 128 # bins per frequency band
    bands = N//(2*bins) # number of frequency bands

    ff, tds_coherence, tds_coefs = noisesub.mcoherence(tar, wit, FS,
                                                       nperseg=2*bands)
                                                    #    noverlap=0, window="boxcar")
    tds_coefs = tds_coefs[0]

    plt.figure()
    plt.suptitle('Coherence of estimated noise with input signal')
    plt.semilogx(ff, tds_coherence, 'b', label='coherence')
    plt.semilogx([10, FS/2], [0.05, 0.05], 'r--', label="cutoff")
    plt.xlim([10, FS/2])
    plt.legend(loc='upper left')
    plt.grid('on')
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Coherence')
    plt.savefig('coherence' + file_end + '.pdf')

    if(signif):
        keep = np.greater(abs(tds_coherence), 0.15)
        tds_coefs = tds_coefs[keep]
        ff_td = ff[keep]
        print 'set', keep.size-np.sum(keep),'coefs to 0 (of', keep.size,').'

    # get the predicted filter frequency response using signal.freqz:
    filt_coefs = coupling_filter()
    Nc = keep.size
    filt_resp = np.ones(Nc)
    for coef in filt_coefs:
        b,a = coef
        w,r = sig.freqz(b,a,worN=Nc)
        filt_resp = filt_resp*np.abs(r)
    freqf = (FS * 0.5 / np.pi) * w
    # We "double pass" the filtering using filtfilt, so we square the filter response
    filt_resp = filt_resp**2
    const = np.mean(np.absolute(tds_coefs)/filt_resp[keep])
    plt.figure()
    plt.suptitle('Transfer function strength')
    # plt.loglog(ff_fd, np.absolute(fds_coefs),'b', label='fd segmented')
    try:
        plt.loglog(ff_td, np.absolute(tds_coefs),'r', label='estimated tf magnitude')
        plt.loglog(freqf, const*filt_resp, 'k', label='actual tf')
        plt.xlim([10,1000])
        plt.ylim(const/4, const*4)
        plt.legend()
        plt.xlabel('Freq (Hz)')
        plt.ylabel('Transfer function magnitude')
        plt.savefig('coefs' + file_end + '.pdf')
    except ValueError:
        print "No significant coefs to plot"

def plot_timeseries_fit(f_vals, f_guess_vals, file_end):
    '''
    Generate plot showing a small segment of the timeseries for the estimated
    and actual added noise. Only useful if f(y1,y2) is not dominated by low
    frequencies (< ~40 Hz).
    '''
    # calculate best fit constant to multiply with f_guess (can be negative)
    # tolerance may be too small for badly approximated functions
    const = scipy.optimize.minimize(
        (lambda c: np.sum((f_vals - c[0]*f_guess_vals-c[1])**2)/f_vals.size),
        [1.0, 0.0],
        method="SLSQP",
        options={'ftol':1e-4*np.std(f_vals)**2})

    if(const.success == False):
        print "Optimization failed:",const.message
    times = np.linspace(0, SEC, FS*SEC, endpoint=False)
    plt.figure()
    plt.plot(times[:FS//10], f_vals[:FS//10], 'b', label = 'f(y1,y2)')
    plt.plot(times[:FS//10], const.x[0]*f_guess_vals[:FS//10]+const.x[1], 'r--', label = 'gp fit')
    plt.legend(fontsize=12)
    plt.xlabel('time (s)')
    plt.xlim([0, 0.1])
    plt.savefig('timeseries_fit' + file_end + '.pdf')

def plot_sub_results(data, tar_orig, wit_guess, f_vals, file_end="", td_seg=True, filt=True):
    '''
    Generate and save plot showing PSDs for background, added noise, and
    subtraction residuals.
    '''
    fs=FS
    y_guess, coh, fit = coherence_sub(tar_orig, wit_guess, td_seg=td_seg)
    y_guess = np.fft.irfft(y_guess, len(data))

    y_guess2, coh2, fit2 = coherence_sub(tar_orig, NOISE_LEVEL*f_vals, td_seg=td_seg)
    y_guess2 = np.fft.irfft(y_guess2, len(data))

    resid = y_guess - data
    resid2 = y_guess2 - data
    nyq = fs/2.0

    Pxx_data, freqs = mlab.psd(data, Fs = fs, NFFT = fs/2, noverlap=fs/4)
    asd_data = np.sqrt(Pxx_data)

    Pxx_tar, freqs = mlab.psd(tar_orig, Fs = fs, NFFT = fs/2, noverlap=fs/4)
    asd_tar = np.sqrt(Pxx_tar)

    Pxx_f, f = mlab.psd(apply_filt(f_vals, filt), Fs = fs, NFFT = fs/2, noverlap=fs/4)
    #freqs2, Pxx_f = sig.welch(apply_filt(f_vals), fs=fs, nperseg=fs/2)
    asd_f = np.sqrt(Pxx_f)

    Pxx_r, freqs = mlab.psd(y_guess, Fs = fs, NFFT = fs/2, noverlap=fs/4)
    asd_r = np.sqrt(Pxx_r)

    Pxx_resid, f = mlab.psd(resid, Fs = fs, NFFT = fs/2, noverlap=fs/4)
    #freqs2, Pxx_resid = sig.welch(resid, fs=fs, nperseg=fs/2)
    asd_resid = np.sqrt(Pxx_resid)

    Pxx_perfect, f = mlab.psd(resid2, Fs = fs, NFFT = fs/2, noverlap=fs/4)
    asd_perfect = np.sqrt(Pxx_perfect)

    worse = np.greater(Pxx_resid, 1.1*Pxx_f)
    print "how many bins worse:", np.sum(worse)
    if(np.sum(worse) <= 20):
        print "at frequencies:",freqs[worse]


    plt.figure()
    plt.loglog(freqs, asd_data, label='background', linewidth=3.5)
    plt.loglog(freqs, asd_f, label='nonlin noise', linewidth=2.5)
    plt.loglog(freqs, asd_tar, label='bg + noise', linewidth=2.5)
    plt.loglog(freqs, asd_r,label='result', linewidth=2.5)
    plt.loglog(freqs, asd_resid,label='residual', linewidth=2.5)
    plt.loglog(freqs, asd_perfect, label='ideal resid', linewidth=2.5)

    plt.xlim([10, fs/2])
    plt.ylim([1e-25, 1e-19])
    plt.grid('on')
    plt.ylabel('ASD (strain/rtHz)')
    plt.legend(loc='upper right', ncol=2)
    plt.xlabel('Freq (Hz)')
    plt.savefig('sub_results' + file_end + '.pdf')

def plot_results(tar, y1, y2, f_guess_vals, data, file_end="", plot_all=True, filt=True):
    '''
    Make several plots to indicate how well the guessed function allows
    subtraction of the added noise, including PSDs, coherence, timeseries fit
    and test data PSD. If plot_all = False, make only the PSD plot with
    training data.
    '''
    # calculate subtraction results and other relevant metrics
    f_vals = f(y1,y2)
    tar_orig = apply_filt(f(y1,y2), filt) + data
    wit_guess = NOISE_LEVEL*f_guess_vals

    # set styles
    plt.style.use('seaborn-colorblind')
    plt.style.use('seaborn-whitegrid')
    plt.style.use('seaborn-talk')

    # main results plot
    plot_sub_results(data, tar, wit_guess, f_vals, file_end, filt=filt)

    if(plot_all):
        coherence_plots(tar, wit_guess, signif=True, file_end=file_end)
        plot_timeseries_fit(f_vals, f_guess_vals, file_end)
