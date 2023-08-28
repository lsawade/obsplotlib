import typing as tp
import obspy
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel


def add_traveltime(stream: obspy.Stream, event_depth_in_m: float = 0.0,
                   origin_time: obspy.UTCDateTime | None = None,
                   phase: str = 'P', vlove=4.4, vrayleigh=3.7, orbit=1,
                   model: str = "ak135", return_filtered: bool = False):

    # Sort the stream
    newstream = stream.copy()
    newstream.sort(keys=['distance', 'network', 'station'])

    # For filtering of the stream if desired. Only relevant for
    # body waves and if return_filtered is True.
    poppable = []

    if phase.lower() == 'love':

        for _i, _tr in enumerate(stream):
            dist = (orbit - 1) * 180.0 + ((orbit-1) % 2) * \
                180 + (-1)**(orbit-1) * _tr.stats.distance
            distkm = dist * 111.11  # km/deg
            _tr.stats.traveltime = distkm/vlove
            _tr.stats.tt_correction_type = 'linear'
            _tr.stats.tt_correction_velocity = vlove
            _tr.stats.label = f'L{orbit}({vlove:.2f} km/s)'

            if origin_time is not None:
                _tr.stats.origin_time = origin_time

    elif phase.lower() == 'rayleigh':

        for _i, _tr in enumerate(stream):
            dist = (orbit - 1) * 180.0 + ((orbit-1) % 2) * \
                180 + (-1)**(orbit-1) * _tr.stats.distance

            distkm = dist * 111.11  # km/deg
            _tr.stats.traveltime = distkm/vrayleigh
            _tr.stats.tt_correction_type = 'linear'
            _tr.stats.tt_correction_velocity = vrayleigh
            _tr.stats.label = f'R{orbit}({vrayleigh:.2f} km/s)'

            if origin_time is not None:
                _tr.stats.origin_time = origin_time

    else:

        # Initialize Taup model
        model = TauPyModel(model=model)

        for _i, _tr in enumerate(stream):

            arrivals = model.get_travel_times(
                source_depth_in_km=event_depth_in_m/1000,
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

                if origin_time is not None:
                    _tr.stats.origin_time = origin_time

    # Remove traces without traveltime
    if return_filtered:
        newstream = stream.copy()

        # Remove traces without traveltime
        for _i in poppable[::-1]:
            newstream.pop(_i)

        return newstream
