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
import numpy as np
import pylab as plt

def main():
    file_pattern = "/data2/dfm/lucky/skemer/*.fits"
    entries = glob.glob(file_pattern)
    image_list = [os.path.abspath(e) for e in sorted(entries)]
    assert len(self.image_list) > 0, \
        "There are no files matching '{0}'".format(file_pattern)
    print image_list
    return None

if __name__ == '__main__':
    main()
