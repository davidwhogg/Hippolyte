'''
This file is part of the Hippolyte project.
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
    return the `(nx, ny)` mean image, the `(K, nx, ny)` eigenimages,
    all the `(N, K)` coefficients, and all `(N, nx, ny)` reconstructed
    images (from the `K` eigenimages).
    """
    N, nx, ny = images.shape
    meanimage = np.mean(images, axis=0)
    datamatrix = (images - meanimage[None, :, :]).reshape(N, nx * ny)
    u, s, vh = la.svd(datamatrix, full_matrices=False)
    # check component ordering
    assert np.all(s[:-1] - s[1:] >= 0.)
    if plotname is not None:
        plt.clf()
        plt.plot(np.arange(len(s)-1) + 1., s[:-1], 'k-')
        plt.xlabel('eigenimage')
        plt.ylabel('variance contribution')
        plt.loglog()
        plt.savefig(plotname)
    eigenimages = vh[:K, :].reshape(K, nx, ny)
    # not completely sure of the "transpose" of the following line
    coefficients = u[:, :K] * s[None, :K]
    reconstructs = np.dot(coefficients, vh[:K, :]).reshape(N, nx, ny) + meanimage[None, :, :]
    return meanimage, eigenimages, coefficients, reconstructs

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
    print "main: Starting PCA..."
    meanimage, eigenimages, coefficients, recons = pca(images, 16, plotname="pca.png")
    print "main: ...done PCA"
    print meanimage.shape, eigenimages.shape, coefficients.shape, recons.shape
    v1, v2 = np.percentile(images, [0.01, 99])
    plt.figure(figsize=(8,4))
    plt.gray()
    for j in [0, 1, 2, 3, 5, 7, 11, 13, 17]:
        plt.clf()
        plt.subplots_adjust(left=0.05, right=0.95, wspace=0.05)
        plt.subplot(1,2,1)
        imj = images[j, 32:32+64, 32:32+64]
        rej = recons[j, 32:32+64, 32:32+64]
        plt.imshow(imj, vmin=v1, vmax=v2, interpolation="nearest")
        plt.subplot(1,2,2)
        plt.imshow(rej, vmin=v1, vmax=v2, interpolation="nearest")
        plt.savefig("recon-%03d.png" % j)
    return None

if __name__ == '__main__':
    main()
