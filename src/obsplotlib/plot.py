import matplotlib as mpl

from .attribdict import AttribDict
from .process import process
from .section import section, section_multiple_comp
from .seismogram import trace, station
from .stream_utils import param_in_streams, stream_max, stream_min, \
    attach_geometry, copy_geometry, select_intersection
from .traveltime import add_traveltime
from .utils import plot_label, add_header
from .window import Window

__all__ = ['AttribDict', 'section', 'process', 'plot_label',
           'param_in_streams', 'stream_max', 'stream_min', 'process',
           'trace', 'station', 'attach_geometry', 'copy_geometry',
           'add_traveltime', 'section_multiple_comp', 'add_header',
           'select_intersection', 'Window']

mpl.rcParams["font.family"] = "monospace"
