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
    ax.set_yticklabels("")
    ax.set_xticklabels("")
    #    asinh = lambda img, mu, sigma, f: f * np.arcsinh((img - mu) / sigma) + 0.2
    #    plotimage = asinh(image[64-32:64+32,64-32:64+32], np.median(image), sigma, foo)
    plotimage = (image[64-32:64+32,64-32:64+32] - np.median(image)) / sigma
    return ax.imshow(plotimage, vmin=-40., vmax=200.0, interpolation="nearest")

def hogg_savefig(fn):
    print "saving {0}".format(fn)
    return plt.savefig(fn)

def main():
    file_pattern = "/data2/dfm/lucky/skemer/*.fits"
    entries = glob.glob(file_pattern)
    files = [os.path.abspath(e) for e in sorted(entries)]
    assert len(files) > 0, \
        "There are no files matching '{0}'".format(file_pattern)
    nn = 6
    plt.gray()
    figsize = (10,10)
    plt.figure(figsize=figsize)
    for i,fn in enumerate(files):
        hdulist = pyfits.open(fn)
        image = hdulist[0].data
        hdulist.close()
        print i, fn, image.shape
        if i == 0:
            sigma = estimate_sigma(image)
        plt.subplot(nn, nn, (i % (nn * nn)) + 1)
        hogg_imshow(plt.gca(), image, sigma, 0.1)
        if (i + 1) % (nn * nn) == 0:
            hogg_savefig("{0}.png".format(os.path.basename(fn)))
            plt.figure(figsize=figsize)
            assert False
    return None

if __name__ == '__main__':
    main()
