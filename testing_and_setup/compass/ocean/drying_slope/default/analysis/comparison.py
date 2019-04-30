#!/usr/bin/env python
"""

Drying slope comparison betewen MPAS-O, analytical, and ROMS results from

Warner, J. C., Defne, Z., Haas, K., & Arango, H. G. (2013). A wetting and
drying scheme for ROMS. Computers & geosciences, 58, 54-61.

Phillip J. Wolfram and Zhendong Cao
04/30/2019

"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# render statically by default
plt.switch_backend('agg')

def setup_fig():
    fig, ax = plt.subplots(nrows=2,ncols=1, sharex=True, sharey=True)
    fig.text(0.04, 0.5, 'Channel depth (m)', va='center', rotation='vertical')
    fig.text(0.5, 0.02, 'Along channel distance (km)', ha='center')

def setup_subplot():
    plt.xlim(0,25)
    plt.ylim(-1, 11)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    x = np.linspace(0,25,100)
    y = 10.0/25.0*x
    plt.plot(x, y, 'k-', lw=3)

def upper_plot():
    plt.subplot(2,1,1)
    plt.gca().set_xticklabels([])
    setup_subplot()

def lower_plot():
    plt.subplot(2,1,2)
    setup_subplot()


setup_fig()

############### subplot r = 0.0025 ###############
upper_plot()



############### subplot r = 0.01   ###############
lower_plot()

# data from MPAS-O on boundary
#ds = xr.open_mfdataset('output.nc')
#ds.ssh.where(ds.tidalInputMask).mean('nCells').plot(marker='o', label='MPAS-O')

#plt.legend()

plt.suptitle('Drying slope comparison between MPAS-O, ROMS, and analytical')

plt.savefig('dryingslopecomparison.png')
