
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
        event_origin_time: obspy.UTCDateTime | None = None,
        lw: list | float = 1.0,
        ls: list | str = '-',
        absmax: float | None = None,
        normalization_type: str | int = 'all',
        plot_labels: bool = True,
        legend: bool = True,
        **kwargs):

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

    # Set some label arguments and time normalization if event_origin_time
    # is given
    if event_origin_time is not None:

        # Get total time delta either from the limits or from the endtimes
        # of the traces compared to the event origin time
        if limits is not None:
            delta = limits[1] - limits[0]
        else:
            latest_endtime: obspy.UTCDateTime = max(
                [tr.stats.endtime for tr in traces])
            delta = latest_endtime - event_origin_time

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
        # If no event_origin_time is given, use matplotlib time
        xlabel = 'Time'
        xdiv = 1

    # Given the divisor, get the corrected x axis limits
    if limits is not None and event_origin_time is not None:
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
    if event_origin_time is not None:
        time_args = dict(type='relative', reftime=event_origin_time)
    else:
        time_args = dict(type='matplotlib')

    # Loop over set of traces
    for _j, _tr in enumerate(traces):

        # Get time vector
        t = _tr.times(**time_args)

        # Normalize time if event_origin_time is given
        if event_origin_time is not None:
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
                   location=3, fontsize='medium', box=False)

        # Plot label if absmax automatically determined
        if not absmax_given:
            plot_label(
                ax, f'|A|max: {absmax:.5g} m', dist=0,
                fontsize='xx-small', box=False)

    # Set limits
    if limits is not None:

        # Set datetime limits if no event origin time is given
        if event_origin_time is None:
            ax.set_xlim([lim.datetime for lim in limits])
        else:
            ax.set_xlim(limits)

    # Set labels and formatting there of
    if event_origin_time is None:

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
                   ncol=3, fontsize='x-small')

    return ax


def station(streams: tp.List[obspy.Stream] | obspy.Stream, *args,
            components: str = "ZRT",
            headerdict: dict | None = None,
            **kwargs):

    # Make sure a list of streams is parsed
    if isinstance(streams, obspy.Stream):
        streams = [streams]

    fig = plt.figure()
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

        # Format x axis to have the date
        if _i < len(components) - 1:
            ax.spines.bottom.set_visible(False)
            ax.tick_params(bottom=False)

    plt.subplots_adjust(hspace=0.0)
