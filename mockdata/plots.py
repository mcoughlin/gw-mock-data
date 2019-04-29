import os
import numpy as np
import scipy.signal as sig
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

__all__ = ['plot_asd', 'plot_coherence']

fs_d = 2048
tfft_d = 8


def _save_incremented(fig, prefix):
    ii = 0
    while os.path.exists('{}_{}.pdf'.format(prefix, ii)):
        ii += 1
    return fig.savefig('{}_{}.pdf'.format(prefix, ii))


def plot_asd(bg, tar, est=None, fs=fs_d, tfft=tfft_d, save=True, model=None):

    if model is None:
        title = 'Mock Data ASD Comparison'
        outfile = 'asd_comparison.pdf'
    else:
        title = 'ASDs of {} Mock Data'.format(model)
        outfile = '{}_asd_comparison.pdf'.format(model)

    nperseg = fs * tfft

    ff, pbg = sig.welch(bg, fs=fs, nperseg=nperseg)
    _, ptar = sig.welch(tar, fs=fs, nperseg=nperseg)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.loglog(ff, np.sqrt(ptar), label='Subtraction Target', rasterized=True)
    ax.loglog(ff, np.sqrt(pbg), label='Background', rasterized=True)

    if est is not None:
        sub = tar - est
        mismatch = bg - sub
        _, psub = sig.welch(sub, fs=fs, nperseg=nperseg)
        _, pmis = sig.welch(mismatch, fs=fs, nperseg=nperseg)
        ax.loglog(ff, np.sqrt(psub), label='Ideal Residual', rasterized=True)
        ax.loglog(ff, np.sqrt(pmis), label='Mismatch', rasterized=True)

    ax.grid(True, which='major')
    ax.grid(True, which='minor', alpha=0.3)
    ax.set_xlim([3, ff[-1]])
    ax.set_ylim([np.min(np.sqrt(pbg)) / 3, np.max(np.sqrt(ptar)) * 3])
    ax.set_xlabel('Freq [Hz]')
    ax.set_ylabel('ASD [1/rtHz]')
    ax.set_title(title)
    ax.legend(frameon=True, framealpha=0.7, loc='best')

    print('Plotting to {}'.format(outfile))
    plt.savefig(outfile, dpi=200)

    return fig, ax


def plot_coherence(tar,
                   est,
                   fs=fs_d,
                   tfft=tfft_d,
                   logy=False,
                   save=False,
                   title=None):
    if title is None:
        title = 'Coherence of Target with Estimate'

    nperseg = fs * tfft

    ff, C = sig.coherence(tar, est, fs=fs, nperseg=nperseg)

    fig, ax = plt.subplots()

    ax.plot(ff, C, label='Coherence with Estimate')

    ax.set_xscale('log')
    ax.set_xlim([10, fs // 2])
    if logy:
        ax.set_yscale('log')
        ax.set_ylim([1e-3, 1.05])
    else:
        ax.set_ylim([0, 1.05])

    ax.set_xlabel('Freq [Hz]')
    ax.set_ylabel('Coherence')
    ax.set_title(title)

    if save:
        _save_incremented(fig, 'coherence')
    else:
        plt.draw_if_interactive()

    return fig, ax
