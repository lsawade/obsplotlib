import matplotlib as mpl


from .attribdict import (
    AttribDict)
from .process import (
    process)
from .section import (
    section,
    section_multiple_comp)
from .seismogram import (
    station,
    trace)
from .signal import (
    L2,
    X,
    Xratio)
from .stream_utils import (
    attach_geometry,
    copy_geometry,
    copy_trace_param,
    make_measurements,
    param_in_streams,
    select_intersection,
    stream_max,
    stream_min)
from .traveltime import (
    add_traveltime)
from .utils import (
    add_header,
    plot_label)
from .window import (
    Window)

__all__ = ['AttribDict',
           'process',
           'section',
           'section_multiple_comp',
           'station',
           'trace',
           'L2',
           'X',
           'Xratio',
           'attach_geometry',
           'copy_geometry',
           'copy_trace_param',
           'make_measurements',
           'param_in_streams',
           'select_intersection',
           'stream_max',
           'stream_min',
           'add_traveltime',
           'add_header',
           'plot_label',
           'Window']

mpl.rcParams["font.family"] = "monospace"
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'monospace'
mpl.rcParams['mathtext.bf'] = 'monospace:bold'
