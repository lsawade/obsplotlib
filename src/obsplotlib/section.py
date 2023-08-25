import numpy as np
import matplotlib.dates as mdates
import matplotlib.axes
import matplotlib.pyplot as plt
import typing as tp
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
import obspy
from .utils import plot_label

from . import stream_utils as su


def section(streams: tp.List[obspy.Stream], *args,
            event_origin_time: obspy.UTCDateTime | None = None,
            ax: matplotlib.axes.Axes | None = None, comp='Z',
            limits: tp.Tuple[obspy.UTCDateTime] | None = None,
            scale: float = 1.0,
            colors=['k', 'r', 'b'],
            labels=['Observed', 'Synthetic', 'New'],
            legendargs: dict | None = None,
            align: bool = False,
            **kwargs):

    # If axes is None generate new axes
    if ax is None:
        ax = plt.gca()

    # If single stream is given, make it a list
    if isinstance(streams, obspy.Stream):
        streams = [streams]

    # Get a single component
    pstreams = [stream.select(component=comp).copy() for stream in streams]

    # Number of datasets
    Nstreams = len(pstreams)

    # Make more colors if needed and print warning
    if len(colors) < Nstreams:
        cmap = plt.get_cmap('rainbow')
        colors = ['k', *cmap(np.linspace(0, 1, Nstreams-1, endpoint=True))]

    # Number of traces
    Ntraces = len(pstreams[0])

    # Sort streams
    for _st in pstreams:
        _st.sort(keys=['distance', 'network', 'station'])

    # If align is True, check all traces for traveltime and window
    if align:

        if event_origin_time is None:
            raise ValueError('event_origin_time must be given if align=True \n'
                             'since traveltime is with respect to origin')

        if not su.param_in_streams(streams, 'traveltime', dtype=float):
            raise ValueError('traveltime not in streams')

    # Get scaling
    if limits is not None:
        slicestart = limits[0]
        sliceend = limits[1]
        maxs = []
        for _st in pstreams:
            maxs.append(np.max([np.max(np.abs(_tr.copy().slice(slicestart, sliceend).data))
                                for _tr in _st]))
        absmax = np.max(maxs)
    else:
        maxs = []
        for _st in pstreams:
            maxs.append(np.max([np.max(np.abs(_tr.data)) for _tr in _st]))
        absmax = np.max(maxs)

    # Plot overall max amplitude label
    plot_label(ax, f'max|u|: {absmax:.5g} m',
               fontsize='small', box=False, dist=0.0, location=4)

    # Plot component label
    plot_label(ax, f'{comp}', fontweight='bold',
               fontsize='medium', box=False, dist=0.0, location=1)

    # Number of stations
    y = np.arange(1, 1*len(pstreams[0])+1, 1)

    # Set ylabels on the left axis to station info
    ax.set_yticks(
        y, [f"{tr.stats.network}.{tr.stats.station}\n"
            f"D:{tr.stats.distance:>6.2f}\n"
            f"Az: {tr.stats.azimuth:>5.1f}" for tr in pstreams[0]],
        verticalalignment='center',
        horizontalalignment='right',
        fontsize='small')

    # Set y values on the right axis to max values
    ax2 = ax.secondary_yaxis("right")

    # Make ticks contain the absolute max value of each trace
    ticks = []
    for _i in range(Ntraces):
        label = ""
        for _j in range(Nstreams):
            label += f"{labels[_j][0].upper()}:{np.max(np.abs(pstreams[_j][_i].data)):>10.4g}"
            # Add newline if not last trace
            if Nstreams > 1 and _j != Nstreams - 1:
                label += "\n"
        ticks.append(label)

    # Set tick labels
    ax2.set_yticks(y, ticks,
                   verticalalignment='center',
                   horizontalalignment='left',
                   fontsize='small')

    # Remove spine
    ax2.spines.right.set_visible(False)

    # Remove ticks
    ax2.tick_params(left=False, right=False)

    # Set time arguments for the plotting
    if event_origin_time is not None:
        time_args = dict(type='relative', reftime=event_origin_time)
    else:
        time_args = dict(type='matplotlib')

    # Normalize
    for _i in range(Ntraces):

        # Get max values
        absmax = np.max([np.max(np.abs(pstreams[_j][_i].data))
                         for _j in range(Nstreams)])

        for _j in range(Nstreams):

            if _i == 0:
                label = labels[_j]
            else:
                label = None

            plt.plot(
                pstreams[_j][_i].times(**time_args),
                pstreams[_j][_i].data /
                absmax * scale + y[_i], '-',
                *args, c=colors[_j], label=label, **kwargs)

    # Remove all spines
    ax.spines.top.set_visible(False)
    ax.spines.left.set_visible(False)
    ax.spines.right.set_visible(False)
    ax.tick_params(left=False, right=False)

    # Plot legend if wanted by setting the legendargs
    if legendargs is not None:
        plt.legend(**legendargs)

    # Format x axis to have the date
    if event_origin_time is None:
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(
            mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))

        if limits is not None:
            ax.set_xlim([lim.datetime for lim in limits])

        plt.xlabel('Time')
    else:
        ax.set_xlim(limits)
        plt.xlabel('Time since origin (s)')


def section_multiple_comp(
    *args, figsize=(12, 4.5),
    components: str = "NEZ",
    legendargs: dict | None = None,
        **kwargs):

    fig = plt.figure(figsize=figsize)
    axes = []
    for _i, comp in enumerate(components):
        if _i == len(components) - 1:
            print(_i, 'legendargs')
            _legendargs = legendargs
        else:
            _legendargs = None

        ax = plt.subplot(1, len(components), _i+1)

        section(*args, ax=ax, comp=comp, legendargs=_legendargs, **kwargs)

        if _i > 0:
            ax.tick_params(labelleft=False)

        axes.append(ax)

    plt.subplots_adjust(left=0.075, right=0.9, top=0.875, wspace=0.45)

    return fig, axes
