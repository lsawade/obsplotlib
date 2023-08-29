
import obspy
import typing as tp
import numpy as np
import matplotlib.axes
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from .utils import plot_label


def trace(
        traces: tp.List[obspy.Trace] | obspy.Trace, *args,
        ax: matplotlib.axes.Axes | None = None,
        limits: tp.Tuple[obspy.UTCDateTime, obspy.UTCDateTime] | None = None,
        nooffset: bool = False,
        colors: list = ['k', 'r', 'b'],
        labels: tp.List[str] = ['Observed', 'Synthetic', 'New Synthetic'],
        origin_time: obspy.UTCDateTime | None = None,
        lw: list | float = 1.0,
        ls: list | str = '-',
        absmax: float | None = None,
        normalization_type: str | int = 'all',
        plot_labels: bool = True,
        legend: bool = True,
        **kwargs):
    """Plot a single or a list of traces.

    Parameters
    ----------
    traces : tp.List[obspy.Trace] | obspy.Trace
        obspy.Trace or list of traces
    ax : matplotlib.axes.Axes | None, optional
        plot into existing axes, by default None
    limits : tp.Tuple[obspy.UTCDateTime, obspy.UTCDateTime] | None, optional
        set axes limits to UTCDateTimes or floats if origin_time is given,
        by default None
    nooffset : bool, optional
        if you are comparing two traces they are automatically offset, set True
        if you don't want that, by default False
    colors : list, optional
        list of colors for traces, has to be same length or longer than
        the number of traces provided, by default ['k', 'r', 'b']
    labels : tp.List[str], optional
        list of colors for labels, has to be same length or longer than
        the number of traces provided,
        by default ['Observed', 'Synthetic', 'New Synthetic']
    origin_time : obspy.UTCDateTime | None, optional
        plot traces with respect to some origin time, by default None
    lw : list | float, optional
        line width(s) can be single float or list thereof. List has to be same
        length or longer than the number of traces provided, by default 1.0
    ls : list | str, optional
        line styles(s) can be single float or list thereof. List has to be same
        length or longer than the number of traces provided, by default '-'
    absmax : float | None, optional
        normalize traces with respect to a specific absmax, by default None
    normalization_type : str | int, optional
        finds absmax to normalize by with respect trace with this index or with
        respect to 'all' traces, by default 'all'
    plot_labels : bool, optional
        Plot trace label, absmax amplitude, by default True
    legend : bool, optional
        plot a legend or not, by default True

    Returns
    -------
    matplotlib.axes.Axes
        returns axes that the traces have been plotted to

    Raises
    ------
    ValueError
        Wrong normalization type input
    ValueError

    """

    if ax is None:
        ax = plt.gca()

    if isinstance(traces, obspy.Trace):
        traces = [traces]

    if isinstance(lw, float) or isinstance(lw, int):
        lw = [lw]*len(traces)

    if isinstance(ls, str):
        ls = [ls]*len(traces)

    # Check if normalization_type is valid
    if isinstance(normalization_type, int):
        if normalization_type >= len(traces):
            raise ValueError('normalization_type must be smaller than the '
                             'number of streams')
    else:
        if not normalization_type in ('all', 'first'):
            raise ValueError(
                'normalization_type must be either "all" or "first"')

    # Get limits and timedelta
    if limits is not None:
        starttime, endtime = limits

    # Set some label arguments and time normalization if origin_time
    # is given
    if origin_time is not None:

        # Get total time delta either from the limits or from the endtimes
        # of the traces compared to the event origin time
        if limits is not None:
            delta = limits[1] - limits[0]
        else:
            latest_endtime: obspy.UTCDateTime = max(
                [tr.stats.endtime for tr in traces])
            delta = latest_endtime - origin_time

        # After computing the delta decide which unit to use
        if delta > 1800:
            xlabel = 'Time since origin [min]'
            xdiv = 60
        elif delta > 10800:
            xlabel = 'Time since origin [h]'
            xdiv = 3600
        else:
            xlabel = 'Time since origin [s]'
            xdiv = 1

    else:
        # If no origin_time is given, use matplotlib time
        xlabel = 'Time'
        xdiv = 1

    # Given the divisor, get the corrected x axis limits
    if limits is not None and origin_time is not None:
        limits = [lim / xdiv for lim in limits]

    # Get max amplitude
    if absmax is None:

        if normalization_type == 'all':

            absmaxs = []
            for tr in traces:

                # Copy the traces
                trc = tr.copy()

                # Only consider amplitude within limits for normalization
                if limits is not None:
                    trc.trim(starttime=starttime, endtime=endtime)

                # Get max amplitudeof each trace
                absmaxs.append(np.max(np.abs(trc.data)))

        elif normalization_type == 'first':

            # Copy the traces
            trc = traces[0].copy()

            # Only consider amplitude within limits for normalization
            if limits is not None:
                trc.trim(starttime=starttime, endtime=endtime)

            # Get max amplitudeof each trace
            absmaxs = [np.max(np.abs(trc.data))]

        # Get absmax of all traces
        absmax = np.max(absmaxs)

        # absmax_given
        absmax_given = False
    else:
        absmax_given = True

    # Define offset for two traces
    if len(traces) != 2 or nooffset:
        absmax_off = 0.0
    else:
        absmax_off = 0.1 * absmax

    # Set time arguments for the plotting
    if origin_time is not None:
        time_args = dict(type='relative', reftime=origin_time)
    else:
        time_args = dict(type='matplotlib')

    # Loop over set of traces
    for _j, _tr in enumerate(traces):

        # Get time vector
        t = _tr.times(**time_args)

        # Normalize time if origin_time is given
        if origin_time is not None:
            t /= xdiv

        # plot trace
        plt.plot(t, _tr.data + absmax_off * (-1)**(_j),
                 ls[_j], *args, c=colors[_j], lw=lw[_j], label=labels[_j],
                 **kwargs)

    # Axis limits and indicator
    ax.set_ylim(-1.15*absmax, 1.15*absmax)

    if plot_labels:
        network = traces[0].stats.network
        station = traces[0].stats.station
        component = traces[0].stats.channel[-1]
        plot_label(ax, f'{network}.{station}.{component}', dist=0.025,
                   location=3, fontsize='small', box=False)

        # Plot label if absmax automatically determined
        if not absmax_given:
            plot_label(
                ax, f'|A|max: {absmax:.5g} m', dist=0.025,
                fontsize='small', box=False)

    # Set limits
    if limits is not None:

        # Set datetime limits if no event origin time is given
        if origin_time is None:
            ax.set_xlim([lim.datetime for lim in limits])
        else:
            ax.set_xlim(limits)

    # Set labels and formatting there of
    if origin_time is None:

        ax.xaxis_date()
        ax.xaxis.set_major_formatter(
            mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))

        # Rotate and align the tick labels so they look better.
        ax.figure.autofmt_xdate()

        # Add xlabel
        plt.xlabel(xlabel)

    else:

        # Add xlabel
        plt.xlabel(xlabel)

    # Unnecessary spines and ticks
    ax.tick_params(labelleft=False, left=False)
    ax.spines.right.set_visible(False)
    ax.spines.left.set_visible(False)
    ax.spines.top.set_visible(False)

    # Add legend
    if legend:
        plt.legend(frameon=False, loc='upper right',
                   ncol=3)

    return ax


def station(streams: tp.List[obspy.Stream] | obspy.Stream, *args,
            components: str = "ZRT",
            **kwargs):
    """Plots given set of components of stream(s). Is a wrapper around
    .trace() so :func:`obsplotlib.seismogram.trace`. Streams should only
    contain traces of a single station. Otherwise, result may be unpredictable.

    Parameters
    ----------
    streams : tp.List[obspy.Stream] | obspy.Stream
        Stream or list of Streams to compare.
    components : str, optional
        which components to plot, by default "ZRT"

    Returns
    -------
    list of matplotlib.axes.Axes
        Each axes contains a component. From top to bottom.

    See also
    --------
    :func:`obsplotlib.seismogram.trace`

    """

    # Make sure a list of streams is parsed
    if isinstance(streams, obspy.Stream):
        streams = [streams]

    axes = []

    for _i, comp in enumerate(components):

        if _i == 0:
            _legend = True
        else:
            _legend = False

        # Get the traces for the component
        traces = [st.select(component=comp)[0] for st in streams]

        # Plot the trace
        ax = plt.subplot(len(components), 1, _i+1)

        trace(traces, *args, ax=ax, legend=_legend, plot_labels=False,
              **kwargs)

        # Plot component label
        plot_label(ax, comp, dist=0.025, fontsize='medium', box=False,
                   location=13)

        # Format x axis to have the date
        if _i < len(components) - 1:
            ax.spines.bottom.set_visible(False)
            ax.tick_params(bottom=False)

        axes.append(ax)

    plt.subplots_adjust(hspace=0.0)

    return axes
