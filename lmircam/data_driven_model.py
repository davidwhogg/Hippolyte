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
import scipy.linalg as la
import pylab as plt

def read_data(files):
    """
    Read images in list `files` and return in the form of one huge
    data block.
    """
    images = None
    for i,fn in enumerate(files):
        hdulist = pyfits.open(fn)
        image = hdulist[0].data
        hdulist.close()
        if images is None:
            images = np.zeros((len(files), image.shape[0], image.shape[1]))
            print images.shape
        images[i,:,:] = image
    return images

def pca(images, K, plotname=None):
    """
    Run PCA on a three-dimensional `(N, nx, ny)` data block `images`,
    return the `(nx, ny)` mean image, the `(K, nx, ny)` eigenimages
    and all the `(N, K)` coefficients.
    """
    N, nx, ny = images.shape
    meanimage = np.mean(images, axis=0)
    print meanimage.shape
    datamatrix = (images - meanimage[None, :, :]).reshape(N, nx * ny)
    u, s, vh = la.svd(datamatrix, full_matrices=False)
    print u.shape, s.shape, vh.shape
    print s[0:10], s[-10:]
    # check component ordering
    assert np.all(s[:-1] - s[1:] >= 0.)
    if plotname is not None:
        plt.clf()
        plt.plot(np.arange(len(s)), s, 'k-')
        plt.xlabel('eigenimage')
        plt.ylabel('variance contribution')
        plt.semilogy()
        plt.savefig(plotname)
    eigenimages = vh[:K, :].reshape(K, nx, ny)
    coefficients = u[:, :K]
    return meanimage, eigenimages, coefficients

def main():
    """
    Do everything.
    """
    file_pattern = "/data2/dfm/lucky/skemer/*.fits"
    entries = glob.glob(file_pattern)
    files = [os.path.abspath(e) for e in sorted(entries)]
    assert len(files) > 0, \
        "There are no files matching '{0}'".format(file_pattern)
    images = read_data(files)
    meanimage, eigenimages, coefficients = pca(images, 4, plotname="pca.png")
    print meanimage.shape, eigenimages.shape, coefficients.shape

if __name__ == '__main__':
    main()
