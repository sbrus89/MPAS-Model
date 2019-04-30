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
import pandas as pd

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


def plot_data(rval='0.0025', dtime='0.05', datatype='analytical', *args, **kwargs):
    datafile = 'data/r' + rval + 'd' + dtime + '-' + datatype + '.csv'
    data = pd.read_csv(datafile, header=None)
    measured=plt.scatter(data[0], data[1], *args, **kwargs)


def plot_dataset(rval):
    times = ['0.50', '0.05', '0.40', '0.15', '0.30', '0.25']
    for ii, dtime in enumerate(times):
        plot_data(rval=rval, dtime = dtime, datatype = 'analytical',
                  marker = '.', color = 'b', label='analytical')
        plot_data(rval=rval, dtime = dtime, datatype = 'roms',
                  marker = '.', color = 'g', label='ROMS')
        if ii == 0:
            plt.legend(frameon=False, loc='lower left')
            place_time_labels(times)
            plt.text(0.5, 5, 'r = ' + str(rval))


def place_time_labels(times):
    locs = [9.3, 7.2, 4.2, 2.2, 1.2, 0.2]
    for atime, ay in zip(times, locs):
        plt.text(25.2, ay, atime + ' days', size=8)


def main():
    setup_fig()

    ############### subplot r = 0.0025 ###############
    upper_plot()
    plot_dataset(rval='0.0025')

    ############### subplot r = 0.01   ###############
    lower_plot()
    plot_dataset(rval='0.01')


    # data from MPAS-O on boundary
    #ds = xr.open_mfdataset('output.nc')
    #ds.ssh.where(ds.tidalInputMask).mean('nCells').plot(marker='o', label='MPAS-O')


    plt.suptitle('Drying slope comparison between MPAS-O, analytical, and ROMS')

    plt.savefig('dryingslopecomparison.png')


if __name__ == "__main__":
    main()
