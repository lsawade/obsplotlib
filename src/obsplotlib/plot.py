import matplotlib as mpl

from .attribdict import AttribDict
from .seismogram import trace, station
from .section import section
from .process import process
from .utils import plot_label
from .stream_utils import param_in_streams, stream_max, stream_min, \
    attach_geometry, copy_geometry

__all__ = ['AttribDict', 'seismogram', 'section', 'process', 'plot_label',
           'param_in_streams', 'stream_max', 'stream_min', 'process',
           'trace', 'station', 'attach_geometry', 'copy_geometry']

mpl.rcParams["font.family"] = "monospace"
