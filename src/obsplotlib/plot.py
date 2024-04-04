import matplotlib as mpl


from .attribdict import AttribDict
from .data import download_data
from .dataclasses import SCARDECSTF, CMTSOLUTION
from .frechet import frechet
from .process import process
from .section import section, section_multiple_comp
from .seismogram import station, trace
from .signal import L2, X, Xratio, xcorr
from .stream_utils import (
    attach_geometry,
    copy_geometry,
    copy_trace_param,
    make_measurements,
    param_in_streams,
    select_intersection,
    stream_max,
    stream_min,
)
from .traveltime import add_traveltime, get_arrivals, plot_arrivals
from .utils import add_header, plot_label, axes_from_axes
from .window import Window

__all__ = [
    "AttribDict",
    "download_data",
    "SCARDECSTF",
    "CMTSOLUTION",
    "frechet",
    "process",
    "section",
    "section_multiple_comp",
    "station",
    "trace",
    "L2",
    "X",
    "xcorr",
    "Xratio",
    "attach_geometry",
    "copy_geometry",
    "copy_trace_param",
    "make_measurements",
    "param_in_streams",
    "select_intersection",
    "stream_max",
    "stream_min",
    "add_traveltime",
    "get_arrivals",
    "plot_arrivals",
    "add_header",
    "plot_label",
    "axes_from_axes",
    "Window",
]

mpl.rcParams["font.family"] = "monospace"
mpl.rcParams["mathtext.fontset"] = "custom"
mpl.rcParams["mathtext.rm"] = "monospace"
mpl.rcParams["mathtext.bf"] = "monospace:bold"
