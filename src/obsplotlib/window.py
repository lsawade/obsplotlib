# Just reusing the AttribDict from ObsPy
from obspy.core.util.attribdict import AttribDict
from collections import OrderedDict
import obspy
import numpy as np


def bold(string):
    N = string.count(' ')
    string = f'$\\bf{{{string}}}$'
    string = N * ' ' + string
    return string


class Window(AttribDict):

    starttime: obspy.UTCDateTime
    endtime: obspy.UTCDateTime
    startidx: int
    endidx: int
    measurements: OrderedDict = OrderedDict()

    def __init__(self, tr: obspy.Trace, *args, **kwargs):
        """Initializes a window object with a trace object reference.
        That way we can always access the trace object from the window
        and make measurements that way."""
        super().__init__(*args, **kwargs)

        if 'startidx' in kwargs:
            self.startidx = kwargs['startidx']
            self.endidx = kwargs['endidx']
        else:
            self.startidx = self.get_index(tr, self.starttime)
            self.endidx = self.get_index(tr, self.endtime)

        # Correct the starttime to the index
        self.starttime = tr.stats.starttime \
            + self.startidx/tr.stats.sampling_rate

        # Correct the endtime to the index
        self.endtime = tr.stats.starttime \
            + self.endidx/tr.stats.sampling_rate

    def __len__(self):
        return self.endidx - self.startidx

    @property
    def duration(self):
        return self.endtime - self.starttime

    def get_index(self, tr: obspy.Trace, time: obspy.UTCDateTime):

        # Get all UTCDateTime objects from the trace
        times = tr.times(type="utcdatetime")

        # Get the index on the traces with closest time to the input time
        return int(np.argmin(np.abs(times-time)))

    def __repr__(self):
        return f'Window(starttime={self.starttime}, endtime={self.endtime})'

    def get_label(self) -> str:

        outstr = ''

        comparison_keys = list(self.measurements.keys())
        measurement_keys = list(self.measurements[comparison_keys[0]].keys())

        N = 4
        outstr += f'{"":{N}}'
        for key in comparison_keys:
            outstr += bold(f'{key:>6s}')

        # ADd new line
        outstr += '\n'

        for mkey in measurement_keys:
            outstr += f'{mkey+":":{N}}'

            for ckey in comparison_keys:
                outstr += f' {self.measurements[ckey][mkey]:>5.2f}'

            outstr += '\n'

        return outstr
