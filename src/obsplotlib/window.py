# Just reusing the AttribDict from ObsPy
from obspy.core.util.attribdict import AttribDict
import obspy
import weakref
import numpy as np


class Window(AttribDict):

    starttime: obspy.UTCDateTime
    endtime: obspy.UTCDateTime
    startidx: int
    endidx: int
    measurements: dict

    def __init__(self, tr: obspy.Trace, *args, **kwargs):
        """Initializes a window object with a trace object reference.
        That way we can always access the trace object from the window
        and make measurements that way."""
        super().__init__(*args, **kwargs)

    def __len__(self):
        return self.endidx - self.startidx

    @property
    def duration(self):
        return self.endtime - self.starttime

    def __repr__(self):
        return f'Window(starttime={self.starttime}, endtime={self.endtime})'

    # def L2(self, other: obspy.Trace) -> float:
    #     """Returns the L2 norm between two Traces."""
    #     if self.trace.stats.sampling_rate != other.stats.sampling_rate:
    #         raise ValueError('Sampling rates must be the same')

    #     if self.trace.stats.starttime != other.stats.starttime:
    #         raise ValueError('starttime must be the same')

    #     if self.trace.stats.endtime != other.stats.endtime:
    #         raise ValueError('endtime must be the same')

    #     return float(0.5 * np.sum(
    #         (self.trace.data[self.startidx:self.endidx]
    #          - other.data[self.startidx:self.endidx])**2))
