'''
This file is part of the hippolyte project.
Copyright 2012 David W. Hogg (NYU) <http://cosmo.nyu.edu/hogg/>.
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':12})
    rc('text', usetex=True)
import glob
import os
import pyfits
import numpy as np
import pylab as plt

def estimate_sigma(scene, nsigma=3.5, tol=0.0):
    img = scene.flatten()
    mask = np.ones(len(img), dtype=bool)
    ms_old = 0.0
    for i in range(500):
        m = np.median(img[mask])
        ms = np.mean((img[mask] - m) ** 2)
        mask = (img - m) ** 2 < nsigma ** 2 * ms
        if i > 1 and np.abs(ms - ms_old) < tol:
            break
        ms_old = ms
    return np.sqrt(ms)

def hogg_imshow(ax, image, sigma, foo):
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_frame_on(False)
    #    asinh = lambda img, mu, sigma, f: f * np.arcsinh((img - mu) / sigma) + 0.2
    #    plotimage = asinh(image[64-32:64+32,64-32:64+32], np.median(image), sigma, foo)
    plotimage = (image[64-32:64+32,64-32:64+32] - np.median(image)) / sigma
    return ax.imshow(plotimage, vmin=-40., vmax=200.0, interpolation="nearest")

def hogg_savefig(fig, fn):
    print "saving {0}".format(fn)
    return fig.savefig(fn)

def get_amplitude_and_phase(image, k = 1.1525):
    ygrid, xgrid = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    cosimage = np.cos(k * xgrid)
    sinimage = np.sin(k * xgrid)
    A = np.sum(image * cosimage)
    B = np.sum(image * sinimage)
    return np.sqrt(A ** 2 + B ** 2), np.arctan2(B, A)

def main():
    file_pattern = "/data2/dfm/lucky/skemer/*.fits"
    entries = glob.glob(file_pattern)
    files = [os.path.abspath(e) for e in sorted(entries)]
    assert len(files) > 0, \
        "There are no files matching '{0}'".format(file_pattern)
    nn = 6
    plt.gray()
    figsize = (10,10)
    fig1 = plt.figure(figsize=figsize)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1, 1, 1)
    for i,fn in enumerate(files):
        hdulist = pyfits.open(fn)
        image = hdulist[0].data
        hdulist.close()
        print i, fn, image.shape
        if i == 0:
            sigma = estimate_sigma(image)
        if i < (nn * nn):
            ax1 = fig1.add_subplot(nn, nn, (i % (nn * nn)) + 1)
            hogg_imshow(ax1, image, sigma, 0.1)
        if (i + 1) == (nn * nn):
            hogg_savefig(fig1, "{0}.png".format(os.path.basename(fn)))
        amp, phase = get_amplitude_and_phase(image)
        if amp < 60000:
            print fn
            for k in 0.01 * np.arange(150):
                a, p = get_amplitude_and_phase(image, k=k)
                print k, a
            assert False
        ax2.plot(phase, amp, 'ko', alpha=0.5)
        if i > 36:
            break
    hogg_savefig(fig2, "amp-phase.png")
    return None

if __name__ == '__main__':
    main()
