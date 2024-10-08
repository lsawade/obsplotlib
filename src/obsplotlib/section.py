import numpy as np
import matplotlib.dates as mdates
import matplotlib.patches as patches
import matplotlib.axes
import matplotlib.pyplot as plt
import typing as tp
import obspy
from .utils import plot_label

from . import stream_utils as su


def section(
    streams: tp.List[obspy.Stream],
    *args,
    origin_time: obspy.UTCDateTime | None = None,
    ax: matplotlib.axes.Axes | None = None,
    comp="Z",
    limits: tp.Tuple[obspy.UTCDateTime] | None = None,
    scale: float = 1.0,
    colors=["k", "r", "b"],
    ls=["-", "-", "-"],
    lw=[0.75, 0.75, 0.75],
    labels=["Observed", "Synthetic", "New"],
    legendargs: dict | None = None,
    align: bool = False,
    absmax: float | None = None,
    window: bool = False,
    plot_geometry: bool = True,
    dist_label: bool = True,
    az_label: bool = True,
    plot_amplitudes: bool = True,
    plot_stations_right: bool = False,
    skip_station: int | None = None,
    invert: bool = False,
    **kwargs,
):
    """Plots a section of seismograms of given stream or stream set. The
    stream(s) need(s) to contain traces which have the distance parameter
    set.

    Parameters
    ----------
    streams : tp.List[obspy.Stream]
        obspy.Stream or list of obspy.Streams
    origin_time : obspy.UTCDateTime | None, optional
        origin time to plot traces agains, by default None
    ax : matplotlib.axes.Axes | None, optional
        axes to plot the section into, by default None
    comp : str, optional
        which component to plot, by default 'Z'
    limits : tp.Tuple[obspy.UTCDateTime] | None, optional
        time limits, by default None
    scale : float, optional
        since the traces have to be normalized, we can scale them after
        normalization, by default 1.0
    colors : list, optional
        list of colors must be equal or longer than list
        of streams, by default ['k', 'r', 'b']
    labels : list, optional
        list of labels must be equal or longer than list
        of streams, by default ['Observed', 'Synthetic', 'New']
    legendargs : dict | None, optional
        a set of arguments to be parsed to ``plt.legend()`` if None, no legend
        is plotted, by default None
    align : bool, optional
        align traces to a traveltime. traces must have stats.traveltime
        parameter and origin time must be given, by default False
    absmax : float | None, optional
        optional value to normalize traces by. If none, automatically
        determined, by default None
    window : bool, optional
        plot windows if available, by default False
    plot_geometry : bool, optional
        plot geometry information, by default True
    plot_stations_right : bool, optional
        plot station information on the right axis, by default False
    plot_amplitudes : bool, optional
        plot amplitudes on the right axis, by default True
    dist_label : bool, optional
        plot distance label, by default True
    az_label : bool, optional
        plot azimuth label, by default True
    plot_amplitudes : bool, optional
        plot amplitudes on the right axis, by default True

    Returns
    -------
    tuple of matplotlib.axes.Axes
        first axes is the main axes that is plotted into which handles
        yticklabels on the right, the second axes handles yticklabels
        on the right.

    Raises
    ------
    ValueError
        if ``align=True``, but no ``origin_time`` given
    ValueError
        if ``align=True``, but the traces' stats objects do not have the
        traveltime attribute.
    """

    # If axes is None generate new axes
    if ax is None:
        ax = plt.gca()

    # If single stream is given, make it a list
    if isinstance(streams, obspy.Stream):
        streams = [streams]

    # Fix linestyle for single value.
    if isinstance(ls, str):
        ls = [ls] * len(streams)

    # Get a single component
    pstreams = [stream.select(component=comp).copy() for stream in streams]

    # Number of datasets
    Nstreams = len(pstreams)

    # Make more colors if needed and print warning
    if len(colors) < Nstreams:
        cmap = plt.get_cmap("rainbow")
        colors = ["k", *cmap(np.linspace(0, 1, Nstreams - 1, endpoint=True))]

    # Number of traces
    Ntraces = len(pstreams[0])

    # Sort streams
    for _st in pstreams:
        _st.sort(keys=["distance", "network", "station"])

    # If align is True, check all traces for traveltime and window
    if align:
        if origin_time is None:
            raise ValueError(
                "origin_time must be given if align=True \n"
                "since traveltime is with respect to origin"
            )

        if not su.param_in_streams(streams, "traveltime", dtype=float):
            raise ValueError("traveltime not in streams")

    # Get scaling
    if limits is not None:
        slicestart = limits[0]
        sliceend = limits[1]
        maxs = []

        # Get abs maxs for all streams in the windows specificed by limits
        for _st in pstreams:
            streammaxs = []

            for _tr in _st:
                if origin_time is not None:
                    slicestart = origin_time + limits[0]
                    sliceend = origin_time + limits[1]

                    if align:
                        tt = _tr.stats.traveltime
                        slicestart += tt
                        sliceend += tt

                else:
                    slicestart = limits[0]
                    sliceend = limits[1]

                streammaxs.append(
                    np.max(np.abs(_tr.copy().slice(slicestart, sliceend).data))
                )

            maxs.append(streammaxs)

        if absmax is None:
            absmax = np.max(maxs)

    else:
        maxs = []
        for _st in pstreams:
            streammaxs = []
            for _tr in _st:
                streammaxs.append(np.max(np.abs(_tr.data)))
            maxs.append(streammaxs)

        if absmax is None:
            absmax = np.max(np.max(maxs))

    # Plot overall max amplitude label
    # plot_label(
    #     ax, f"max|u|: {absmax:.5g} m", fontsize="small", box=False, dist=0.0, location=4
    # )

    # Plot component label
    # plot_label(
    #     ax,
    #     f"{comp}",
    #     fontweight="bold",
    #     fontsize="medium",
    #     box=False,
    #     dist=0.0,
    #     location=1,
    # )

    # Number of stations
    y = np.arange(1, 1 * len(pstreams[0]) + 1, 1)

    if invert:
        y *= -1

    # Define ylabels on the left axis to station info
    ylabels = []
    for _i, tr in enumerate(pstreams[0]):
        if skip_station is not None and np.mod(_i, skip_station) != 0:
            # ylabels.append("\xb7")
            ylabel = ""
        else:
            ylabel = f"{tr.stats.network}.{tr.stats.station}"

            if plot_geometry:
                if hasattr(tr.stats, "distance") and dist_label:
                    ylabel += f"\nD:{tr.stats.distance:>6.2f}"

                if hasattr(tr.stats, "azimuth") and az_label:
                    ylabel += f"\nAz: {tr.stats.azimuth:>5.1f}"

        ylabels.append(ylabel)

    # Set labels
    ax.set_yticks(
        y,
        ylabels,
        verticalalignment="center",
        horizontalalignment="right",
        fontsize="x-small",
    )

    # Set y values on the right axis to max values
    ax2 = ax.secondary_yaxis("right")

    # Remove spine
    ax2.spines.right.set_visible(False)

    # Remove ticks
    ax2.tick_params(left=False, right=False)

    if plot_amplitudes and not plot_stations_right:
        # Make ticks contain the absolute max value of each trace
        ticks = []
        for _i in range(Ntraces):
            label = ""
            for _j in range(Nstreams):
                label += f"{labels[_j][0].upper()}:{maxs[_j][_i]:>10.4g}"
                # Add newline if not last trace
                if Nstreams > 1 and _j != Nstreams - 1:
                    label += "\n"
            ticks.append(label)

        # Set tick labels
        ax2.set_yticks(
            y,
            ticks,
            verticalalignment="center",
            horizontalalignment="left",
            fontsize="x-small",
        )

    elif (skip_station % 2) == 0 and plot_stations_right:
        ticks = []
        for _i, tr in enumerate(pstreams[0]):
            if (
                skip_station is not None
                and np.mod(_i + skip_station / 2, skip_station) != 0
            ):
                # ylabels.append("\xb7")
                ticks.append("")
            else:
                ylabel = f"{tr.stats.network}.{tr.stats.station}"

                if plot_geometry:
                    f"\nD:{tr.stats.distance:>6.2f}"

                    if hasattr(tr.stats, "azimuth"):
                        ylabel += f"\nAz: {tr.stats.azimuth:>5.1f}"

                ticks.append(ylabel)
        # Set tick labels
        ax2.set_yticks(
            y,
            ticks,
            verticalalignment="center",
            horizontalalignment="left",
            fontsize="x-small",
        )

    else:
        # Remove ticks
        ax2.tick_params(left=False, right=False, labelright=False)

    # Set time arguments for the plotting
    if origin_time is not None:
        time_args = dict(type="relative", reftime=origin_time)
    else:
        time_args = dict(type="matplotlib")

    # Factor to add to invert the traces
    if invert:
        invert_factor = 1
    else:
        invert_factor = 1

    # Normalize
    for _i in range(Ntraces):
        for _j in range(Nstreams):
            if _i == 0:
                label = labels[_j]
            else:
                label = None

            if align:
                tt = pstreams[_j][_i].stats.traveltime
                time_args["reftime"] = origin_time + tt

            plt.plot(
                pstreams[_j][_i].times(**time_args),
                pstreams[_j][_i].data / absmax * scale * invert_factor + y[_i],
                *args,
                c=colors[_j],
                ls=ls[_j],
                lw=lw[_j],
                label=label,
                **kwargs,
            )

            # Plot windows if available
            if window and _j == 0:
                # if windowkwargs is None:
                #     windowkwargs = dict(color=colors[_j], alpha=0.5)

                for window in pstreams[_j][_i].stats.windows:
                    if origin_time is not None:
                        windowstart = window.starttime - origin_time
                        windowend = window.endtime - origin_time
                        duration = windowend - windowstart
                    else:
                        windowstart = window.starttime.matplotlib_date
                        windowend = window.endtime.matplotlib_date
                        duration = windowend - windowstart

                    # Vertical window extension
                    windowy0 = y[_i] - 0.6 * scale * invert_factor
                    windowy1 = y[_i] + 0.6 * scale * invert_factor
                    windowdy = windowy1 - windowy0

                    ax.add_patch(
                        patches.Rectangle(
                            (windowstart, windowy0),
                            duration,
                            windowdy,
                            edgecolor="none",
                            zorder=-1,
                            facecolor=[0.9, 0.9, 0.9],
                            clip_on=True,
                        )
                    )

    # Remove all spines
    ax.spines.top.set_visible(False)
    ax.spines.left.set_visible(False)
    ax.spines.right.set_visible(False)
    ax.tick_params(left=False, right=False)

    # Plot legend if wanted by setting the legendargs
    if legendargs is not None:
        plt.legend(**legendargs)

    # Format x axis to have the date
    if origin_time is None:
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(
            mdates.ConciseDateFormatter(ax.xaxis.get_major_locator())
        )

        if limits is not None:
            ax.set_xlim([lim.datetime for lim in limits])

        plt.xlabel("Time")
    else:
        ax.set_xlim(limits)
        plt.xlabel("Time since origin (s)")

    # Set y limits
    if invert:
        ylim = (
            y[-1] - 1.5,
            y[0] + 1.5,
        )
    else:
        ylim = y[0] - 1.5, y[-1] + 1.5
    # ax.set_ylim(ylim)
    # ax2.set_ylim(ylim)

    return ax, ax2


def section_multiple_comp(
    streams: tp.List[obspy.Stream],
    *args,
    components: str = "NEZ",
    legendargs: dict | None = None,
    absmax: float | None = None,
    align: bool = False,
    limits: tp.Tuple[obspy.UTCDateTime] | tp.Tuple[float] | None = None,
    **kwargs,
):
    """Wrapper around section that plit multiple sections into a single figure.

    Parameters
    ----------
    streams : tp.List[obspy.Stream]
        stream or list of streams to be plotted into a single axes.
    components : str, optional
        _description_, by default "NEZ"
    legendargs : dict | None, optional
        _description_, by default None
    absmax : float | None, optional
        _description_, by default None
    align : bool, optional
        _description_, by default False
    limits : tp.Tuple[obspy.UTCDateTime] | tp.Tuple[float] | None, optional
        _description_, by default None

    Returns
    -------
    _type_
        _description_
    """

    axes = []

    for _i, comp in enumerate(components):
        if _i == len(components) - 1:
            print(_i, "legendargs")
            _legendargs = legendargs
        else:
            _legendargs = None

        ax = plt.subplot(1, len(components), _i + 1)

        # Get overall absmax
        if absmax is None:
            # Check if traveltime is needed
            if align and limits is not None:
                traveltime = True
            else:
                traveltime = False

            # Get absmax
            absmax = su.get_max(streams, traveltime=traveltime, limits=limits)

        ax, ax2 = section(
            streams,
            *args,
            ax=ax,
            comp=comp,
            legendargs=_legendargs,
            limits=limits,
            absmax=absmax,
            align=align,
            **kwargs,
        )

        if _i > 0:
            ax.tick_params(labelleft=False)

        axes.append((ax, ax2))

    plt.subplots_adjust(left=0.075, right=0.9, top=0.875, wspace=0.45)

    return axes
