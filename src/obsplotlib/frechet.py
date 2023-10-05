# External
import os
import sys
import obspy
import typing as tp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

# Internal
from .utils import plot_label


def frechet(rp: obspy.Stream | tp.List[obspy.Stream],
            drp: tp.Dict[str, obspy.Stream] | tp.List[tp.Dict[str, obspy.Stream]],
            network, station, component, limits=[None, None], lw=1.0):

    if not (isinstance(rp, list) and isinstance(drp, list)) and \
            not (isinstance(rp, obspy.Stream) and isinstance(drp, dict)):
        raise ValueError('rp and drp must either be list of streams and list of'
                         'dicts of streams or a stream and a dictionary of streams.')

    if isinstance(rp, list) and isinstance(drp, list):
        if len(rp) != len(drp):
            raise ValueError('Both rp, and drp lists have to have the same '
                             'number of elements')

    if isinstance(rp, obspy.Stream):
        rp = [rp, ]

    if isinstance(drp, dict):
        drp = [drp, ]

    # Number of parameters
    N = 1 + len(drp[0].keys())
    keys = list(drp[0].keys())

    # Number of streams
    Nst = len(rp)

    factor = 1.5
    fig, axes = plt.subplots(N, 1, figsize=(8, 4.5), gridspec_kw={
        'height_ratios': [factor, *((1,)*(N-1))]})

    for i in range(N):

        absmaxs = []
        for j in range(Nst):
            if i == 0:
                key = 'Syn'
                tr = rp[j].select(
                    network=network, station=station, component=component)[0]

            else:
                key = keys[i-1]
                tr = drp[j][key].select(
                    network=network, station=station, component=component)[0]

            plt.sca(axes[i])
            plt.plot(tr.times("matplotlib"), tr.data, 'k', lw=lw)

            # Store trace maximum
            absmaxs.append(np.max(np.abs(tr.data)))

        # Get overall absmax for parameter
        absmax = np.max(absmaxs)

        axes[i].set_ylim(-1.01*absmax*factor, 1.01*absmax*factor)
        plot_label(
            axes[i], f'{key}', location=1, dist=0.0,
            fontsize='small', box=False)

        plot_label(
            axes[i], f'max|A|: {absmax:.5g}', location=2, dist=0.001,
            fontsize='x-small', box=False)

        axes[i].xaxis_date()
        axes[i].xaxis.set_major_formatter(
            mdates.ConciseDateFormatter(axes[i].xaxis.get_major_locator()))
        axes[i].tick_params(labelleft=False, left=False)
        axes[i].spines.right.set_visible(False)
        axes[i].spines.left.set_visible(False)
        axes[i].spines.top.set_visible(False)

        axes[i].set_xlim(limits)

        if i == N-1:
            plt.xlabel('Time')
        else:
            axes[i].spines.bottom.set_visible(False)

    fig.autofmt_xdate()
    plt.subplots_adjust(
        left=0.075, right=0.925, bottom=0.15, top=0.95, hspace=0.0)
