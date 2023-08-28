import obspy
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
import numpy as np
import typing as tp


def get_max(streams: tp.List[obspy.Stream] | obspy.Stream | obspy.Trace,
            *args, **kwargs) -> float:

    if isinstance(streams, obspy.Trace):
        return trace_max(streams, *args, **kwargs)
    elif isinstance(streams, obspy.Stream):
        return stream_max(streams, *args, **kwargs)
    elif isinstance(streams, list) and isinstance(streams[0], obspy.Stream):
        return stream_set_max(streams, *args, **kwargs)
    else:
        raise ValueError('streams must be either a Trace, Stream or a list of '
                         'streams.')


def trace_max(
        tr: obspy.Trace,
        absolute: bool = True,
        window: bool = False,
        traveltime: bool = False,
        limits: None | tp.Tuple[obspy.UTCDateTime, obspy.UTCDateTime] |
        tp.Tuple[float, float] = None):
    """Get trace max value given certain parameters

    Parameters
    ----------
    tr : obspy.Trace
        obspy trace to be analyzed
    absolute : bool, optional
        whether to output absolute max, by default True, by default True
    window : bool, optional
        window, by default False
    limits : None | tp.Tuple[obspy.UTCDateTime, obspy.UTCDateTime], optional
        limits, by default None
    """

    if window:

        raise ValueError("``Window`` not yet implemented.")

    elif limits is not None:

        if isinstance(limits[0], obspy.UTCDateTime):
            starttime = limits[0]
            endtime = limits[1]
        elif traveltime:
            starttime = tr.stats.origin_time + tr.stats.traveltime + limits[0]
            endtime = tr.stats.origin_time + tr.stats.traveltime + limits[1]
        else:
            starttime = tr.stats.origin_time + limits[0]
            endtime = tr.stats.origin_time + limits[1]

        _tr = tr.copy().slice(starttime, endtime)

        if absolute:
            return float(np.max(np.abs(_tr.data)))
        else:
            return float(np.max(_tr.data))

    else:
        if absolute:
            return float(np.max(np.abs(tr.data)))
        else:
            return float(np.max(tr.data))


def stream_max(st: obspy.Stream, *args, **kwargs) -> float:
    """Return the (absolute) maximum absolute value of a stream.

    Parameters
    ----------
    st : obspy.Stream
        Stream to be analyzed
    *args, **kwargs : optional
        pass to trace_max

    Returns
    -------
    float
        absolute maximum value of the all traces in stream
    """

    return np.max([trace_max(tr, *args, **kwargs) for tr in st])


def stream_min(st: obspy.Stream, absolute: bool = True) -> float:
    """Return the (absolute) minimum absolute value of a stream.

    Parameters
    ----------
    st : obspy.Stream
        Stream to be analyzed
    abs : bool, optional
        whether to output absolute min, by default True

    Returns
    -------
    float
        absolute minimum value of the all traces in stream
    """

    if absolute:
        return float(np.min([np.min(np.abs(tr.data)) for tr in st]))
    else:
        return float(np.min([np.min(tr.data) for tr in st]))


def stream_set_max(streams: tp.List[obspy.Stream], *args, **kwargs) -> float:
    """Return the maximum absolute value of a list of streams.

    Parameters
    ----------
    streams : list of obspy.Stream
        Streams to be analyzed
    *args, **kwargs : optional
        pass to stream_max

    Returns
    -------
    float
        maximum absolute value of the all traces in streams
    """

    return np.max([stream_max(st, *args, **kwargs) for st in streams])


def stream_set_min(streams: tp.List[obspy.Stream], absolute=True) -> float:
    """Return the minimum absolute value of a list of streams.

    Parameters
    ----------
    streams : list of obspy.Stream
        Streams to be analyzed
    abs : bool, optional
        whether to output absolute min, by default True

    Returns
    -------
    float
        minimum absolute value of the all traces in streams
    """

    return np.min([stream_min(st, absolute=absolute) for st in streams])


def select_intersection(
        streams: tp.List[obspy.Stream], components: str = 'NEZ'):

    Nstreams = len(streams)
    stations = set()

    for st in streams:
        for tr in st:
            stations.add((tr.stats.network, tr.stats.station))

    # Grabbing the actual pairs
    newstreams = [[] for _ in range(Nstreams)]
    for _net, _sta in stations:
        for _comp in components:
            try:
                subtraces = []
                for st in streams:
                    subtraces.append(st.select(
                        network=_net, station=_sta, component=_comp)[0])

                for i in range(Nstreams):
                    newstreams[i].append(subtraces[i].copy())

            except Exception as e:
                print(e)
                print(f'Cant find {_net}.{_sta}..{_comp}')

    return [obspy.Stream(traces=tracelist) for tracelist in newstreams]


def get_station_params(st: obspy.Stream):
    """Get unique station names, latitudes and longitudes from a stream.

    Parameters
    ----------
    st : Stream
        obspy.Stream with station information in stats

    Returns
    -------
    tuple of lists
        networks, stations, latitudes, longitudes
    """

    # Initialize
    networks, stations, latitudes, longitudes = [], [], [], []

    # Get unique stations
    station_set = set()
    for tr in st:
        station_set.add((tr.stats.network, tr.stats.station))

    # Append to list of networks and stations, latitudes and longitudes
    for network, station in station_set:
        networks.append(network)
        stations.append(station)
        tr = st.select(network=network, station=station)[0]
        latitudes.append(tr.stats.latitude)
        longitudes.append(tr.stats.longitude)

    return networks, stations, latitudes, longitudes


def attach_geometry(st: obspy.Stream,
                    event_latitude: float, event_longitude: float,
                    inv: obspy.Inventory | None = None):
    """Given event coordinates and an inventory, attach station coordinates,
    distance and backzimuth to each trace. If the inventory is not provided,
    the station coordinates are taken from the stats in each Trace of the
    Stream.

    Parameters
    ----------
    st : Stream
        input stream
    event_latitude : float
        event latitude
    event_longitude : float
        event longitude
    inv : Inventory | None, optional
        inventory, by default None
    """

    if inv is not None:
        # Initialize
        networks, stations, latitudes, longitudes = [], [], [], []

        for network in inv:
            for station in network:
                networks.append(network.code)
                stations.append(station.code)
                latitudes.append(station.latitude)
                longitudes.append(station.longitude)
    else:

        # Get unique station name
        networks, stations, latitudes, longitudes = get_station_params(st)

    for network, station, latitude, longitude in \
            zip(networks, stations, latitudes, longitudes):

        try:
            subset = st.select(network=network, station=station)

            latA = event_latitude
            lonA = event_longitude
            latB = latitude
            lonB = longitude

            # Compute distance
            dist_in_m, az_A2B_deg, az_B2A_deg = gps2dist_azimuth(
                latA, lonA, latB, lonB)

            for tr in subset:
                if not hasattr(tr.stats, 'latitude') and \
                        not hasattr(tr.stats, 'longitude'):
                    tr.stats.latitude = latB
                    tr.stats.longitude = lonB
                tr.stats.distance = dist_in_m/1000/(40000/360)
                tr.stats.back_azimuth = az_B2A_deg
                tr.stats.azimuth = az_A2B_deg
        except:
            print(f'traces for station {network}.{station}')


def copy_geometry(st: obspy.Stream, st2: obspy.Stream | tp.List[obspy.Stream]):
    """Copies the station location, distance to event and back azimuth from
    one stream to one other or a list of other streams.

    Parameters
    ----------
    st : Stream
        stream with station information in stats
    st2 : Stream | tp.List[Stream]
        stream to copy station information to
    """

    # Get station names
    networks, stations, _, _ = get_station_params(st)

    # Make sure
    if isinstance(st2, obspy.Stream):
        st2 = [st2,]

    # Loop over stations
    for network, station in zip(networks, stations):

        # Get trace
        tr = st.select(network=network, station=station)[0]

        # Loop over streams
        for stream in st2:

            # Select traces
            sub_stream = stream.select(network=network, station=station)

            for trace in sub_stream:

                trace.stats.latitude = tr.stats.latitude
                trace.stats.longitude = tr.stats.longitude
                trace.stats.distance = tr.stats.distance
                trace.stats.back_azimuth = tr.stats.back_azimuth
                trace.stats.azimuth = tr.stats.azimuth


def copy_trace_param(st: obspy.Stream, st2: obspy.Stream | tp.List[obspy.Stream],
                     param: str):
    """Copies the station location, distance to event and back azimuth from
    one stream to one other or a list of other streams.

    Parameters
    ----------
    st : Stream
        stream with station information in stats
    st2 : Stream | tp.List[Stream]
        stream to a trace parameter to from st

    """

    # Get station names
    networks, stations, _, _ = get_station_params(st)

    # Make sure
    if isinstance(st2, obspy.Stream):
        st2 = [st2,]

    # Loop over stations
    for network, station in zip(networks, stations):

        # Get trace
        stationstream = st.select(network=network, station=station)

        # Loop over streams
        for stream in st2:

            # Select traces
            sub_stream = stream.select(network=network, station=station)

            for trace in sub_stream:

                tr = stationstream.select(component=trace.stats.component)[0]

                setattr(trace.stats, param, getattr(tr.stats, param))


def get_azimuth_distance_traveltime(
        cmt, obs, syn,
        traveltime_window: None | tp.Tuple[str, tp.Tuple[float, float]],
        comp='Z', newsyn=None, vlove=4.4, vrayleigh=3.7, orbit=1,
        model: str = "ak135"):

    # Get a single component
    pobs = obs.select(component=comp).copy()
    psyn = syn.select(component=comp).copy()

    if newsyn is not None:
        pnewsyn = newsyn.select(component=comp).copy()
    else:
        pnewsyn = None

    # Get station event distances, labels
    for _i, (_obs, _syn) in enumerate(zip(pobs, psyn)):

        # Assign lats/lons
        latA = cmt.latitude
        lonA = cmt.longitude
        latB = _syn.stats.latitude
        lonB = _syn.stats.longitude

        # Compute distance
        dist = locations2degrees(latA, lonA, latB, lonB)

        # Compute azimuth
        # dist_in_m, az_A2B_deg, az_B2A_deg = gps2dist_azimuth
        _, az_A2B_deg, az_B2A_deg = gps2dist_azimuth(latA, lonA, latB, lonB)

        # Add info to traces
        _obs.stats.distance = dist
        _syn.stats.distance = dist
        _obs.stats.azimuth = az_A2B_deg
        _syn.stats.azimuth = az_A2B_deg
        _obs.stats.backazimuth = az_B2A_deg
        _syn.stats.backazimuth = az_B2A_deg

        if pnewsyn:
            pnewsyn[_i].stats.distance = dist
            pnewsyn[_i].stats.azimuth = az_A2B_deg
            pnewsyn[_i].stats.backazimuth = az_B2A_deg

    # Sort the stream
    pobs.sort(keys=['distance', 'network', 'station'])
    psyn.sort(keys=['distance', 'network', 'station'])

    if pnewsyn:
        pnewsyn.sort(keys=['distance', 'network', 'station'])

    if traveltime_window is not None:
        # Phase
        phase = traveltime_window[0]
        trims = traveltime_window[1]

        if phase.lower() == 'love':

            for _i, _tr in enumerate(pobs):
                dist = (orbit - 1) * 180.0 + ((orbit-1) % 2) * \
                    180 + (-1)**(orbit-1) * _tr.stats.distance
                distkm = dist * 111.11  # km/deg
                _tr.stats.traveltime = distkm/vlove
                _tr.stats.tt_correction_type = 'linear'
                _tr.stats.tt_correction_velocity = vlove
                _tr.stats.label = f'L{orbit}({vlove:.2f} km/s)'

        elif phase.lower() == 'rayleigh':

            for _i, _tr in enumerate(pobs):
                dist = (orbit - 1) * 180.0 + ((orbit-1) % 2) * \
                    180 + (-1)**(orbit-1) * _tr.stats.distance

                distkm = dist * 111.11  # km/deg
                _tr.stats.traveltime = distkm/vrayleigh
                _tr.stats.tt_correction_type = 'linear'
                _tr.stats.tt_correction_velocity = vrayleigh
                _tr.stats.label = f'R{orbit}({vrayleigh:.2f} km/s)'

        else:

            # Initialize Taup model
            model = TauPyModel(model=model)

            poppable = []
            for _i, _tr in enumerate(pobs):

                arrivals = model.get_travel_times(
                    source_depth_in_km=cmt.depth,
                    distance_in_degree=_tr.stats.distance,
                    phase_list=[phase, ])

                if len(arrivals) == 0:
                    _tr.stats.traveltime = None
                    poppable.append(_i)
                    _tr.stats.label = None
                else:
                    _tr.stats.traveltime = arrivals[0].time
                    _tr.stats.tt_correction_type = 'phase'
                    _tr.stats.tt_correction_velocity = None
                    _tr.stats.label = f'{phase}-Wave'

            for _pop in poppable[::-1]:

                pobs.pop(_pop)
                psyn.pop(_pop)
                if pnewsyn:
                    pnewsyn.pop(_pop)

        for _i, (_obstr, _syntr) in enumerate(zip(pobs, psyn)):

            # Set up trace by trace interpolation
            dt = 0.1
            starttime = cmt.origin_time + _obstr.stats.traveltime + trims[0]
            npts = int(np.round((trims[1] - trims[0])/dt))

            # Interpolation arguments for the inteprolation
            iargs = 1.0/dt,
            ikwargs = dict(method='lanczos', starttime=starttime,
                           npts=npts, time_shift=0.0, a=20)

            # Interpolation
            _obstr.interpolate(*iargs, **ikwargs)
            _syntr.interpolate(*iargs, **ikwargs)

            if pnewsyn:
                pnewsyn[_i].interpolate(*iargs, **ikwargs)

    if newsyn is not None:
        return pobs, psyn, pnewsyn
    else:
        return pobs, psyn


def traveltime_filter_stations(obs1, obs2):
    selection = []
    for _i, (_obs1, _obs1) in enumerate(zip(obs1, obs2)):
        if (_obs1.stats.traveltime is None) and (_obs2.stats.traveltime is None):
            continue
        else:
            selection.append(_i)

    return selection


def param_in_streams(streams: tp.List[obspy.Stream],
                     param: str, dtype=None) -> bool:
    """Check if all traces in a list of streams have a certain parameter in
    the stats object.

    Parameters
    ----------
    streams : list of obspy.Stream
        streams to check
    param : str
        parameter to check for
    dtype: type, optional
        data type of the parameter, by default None
    Returns
    -------
    bool
        True if all traces in all streams have the parameter, False otherwise
    """

    for st in streams:
        for tr in st:
            if not hasattr(tr.stats, param):
                print(f"{tr.id} does not have {param}.")
                return False

            if dtype is not None:
                if not isinstance(getattr(tr.stats, param), dtype):
                    print(f"{tr.id} has {param}, but is not of type {str(dtype)}.")
                    return False

    return True
