import obspy
import typing as tp
import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

def add_header(ax,
               station: str | None = None,
               station_latitude: float | None = None,
               station_longitude: float | None = None,
               station_azimuth: float | None = None,
               station_back_azimuth: float | None = None,
               station_distance_in_degree: float | None = None,
               event: str | None = None,
               event_time: obspy.UTCDateTime | None = None,
               event_latitude: float | None = None,
               event_longitude: float | None = None,
               event_depth_in_km: float | None = None,
               event_Mw: float | None = None,
               bandpass: tp.List[float] | None = None,
               add_newline_event: bool = False,
               add_newline_station: bool = False,
               **kwargs):

    label = ""

    if event is not None:
        label += f"{event}: "

    if add_newline_event:
        label += "\n  "

    if event_time is not None:
        label += f"{event_time.strftime('%Y-%m-%d %H:%M:%S')}  "

    if event_latitude is not None:
        label += f"La {event_latitude:.3f} "

    if event_longitude is not None:
        label += f"Lo {event_longitude:.3f} "

    if event_depth_in_km is not None:
        label += f"Dp {event_depth_in_km:.3f}km "

    if event_Mw is not None:
        label += f"Mw: {event_Mw:.2f}"

    if any([event, event_time, event_latitude, event_longitude,
            event_depth_in_km, event_Mw]):
        label += "\n"

    if station is not None:
        label += f"{station}: "

    if add_newline_station:
        label += "\n  "

    if station_latitude is not None:
        label += f"La {station_latitude:.3f} "

    if station_longitude is not None:
        label += f"Lo: {station_longitude:.3f} "

    if station_azimuth is not None:
        label += f"Az: {station_azimuth:.1f} "

    if station_back_azimuth is not None:
        label += f"Baz: {station_back_azimuth:.1f} "

    if station_distance_in_degree is not None:
        label += f"Dist: {station_distance_in_degree:.3f}dg "

    if bandpass is not None:
        label += f"- Bp: {bandpass[0]:d}-{bandpass[1]:d}s"

    # Set some default kwargs for plot_label
    if 'fontsize' not in kwargs:
        kwargs['fontsize'] = 'medium'

    if 'location' not in kwargs:
        kwargs['location'] = 6

    # Finally add label
    plot_label(ax, label, box=False, **kwargs)


def adjust_spines(ax: Axes, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


def plot_label(ax: Axes, label: str, aspect: float = 1,
               location: int = 1, dist: float = 0.025,
               box: tp.Union[bool, dict] = True, fontdict: dict = {},
               **kwargs):
    """Plots label one of the corners of the plot.
    Plot locations are set as follows::

        17  6  14  7  18
            --------
         5 |1  22  2| 8
        13 |21  0 23| 15
        12 |3  24  4| 9
            --------
        20  11 16 10  19

    Tee dist parameter defines the distance between the axes and the text.

    Parameters
    ----------
    label : str
        label
    aspect : float, optional
        aspect ratio length/height, by default 1.0
    location : int, optional
        corner as described by above code figure, by default 1
    aspect : float, optional
        aspect ratio length/height, by default 0.025
    box : bool
        plots bounding box st. the label is on a background, default true
    Notes
    -----
    :Author:
        Lucas Sawade (lsawade@princeton.edu)
    :Last Modified:
        2021.01.26 18.30
    """
    if type(box) is bool:
        if box:
            boxdict = {'facecolor': 'w', 'edgecolor': 'k'}
        else:
            boxdict = {'facecolor': 'none', 'edgecolor': 'none'}
    else:
        boxdict = box

    # Get aspect of the axes
    aspect = 1.0/get_aspect(ax)

    # Inside
    if location == 0:
        t = ax.text(0.5, 0.5, label,
                horizontalalignment='center', verticalalignment='center_baseline',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 1:
        t = ax.text(dist, 1.0 - dist * aspect, label, horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 2:
        t = ax.text(1.0 - dist, 1.0 - dist * aspect, label,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 3:
        t = ax.text(dist, dist * aspect, label, horizontalalignment='left',
                verticalalignment='bottom', transform=ax.transAxes,
                bbox=boxdict, fontdict=fontdict, **kwargs)
    elif location == 4:
        t = ax.text(1.0 - dist, dist * aspect, label,
                horizontalalignment='right', verticalalignment='bottom',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    # Outside
    elif location == 5:
        t = ax.text(-dist, 1.0, label, horizontalalignment='right',
                verticalalignment='top', transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 6:
        t = ax.text(0, 1.0 + dist * aspect, label, horizontalalignment='left',
                verticalalignment='bottom', transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 7:
        t = ax.text(1.0, 1.0 + dist * aspect, label,
                horizontalalignment='right', verticalalignment='bottom',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 8:
        t = ax.text(1.0 + dist, 1.0, label,
                horizontalalignment='left', verticalalignment='top',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 9:
        t = ax.text(1.0 + dist, 0.0, label,
                horizontalalignment='left', verticalalignment='bottom',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 10:
        t = ax.text(1.0, - dist * aspect, label,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 11:
        t = ax.text(0.0, -dist * aspect, label, horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes,
                bbox=boxdict, fontdict=fontdict, **kwargs)
    elif location == 12:
        t = ax.text(-dist, 0.0, label, horizontalalignment='right',
                verticalalignment='bottom', transform=ax.transAxes,
                bbox=boxdict, fontdict=fontdict, **kwargs)
    elif location == 13:
        t = ax.text(-dist, 0.5, label, horizontalalignment='right',
                verticalalignment='center_baseline', transform=ax.transAxes,
                bbox=boxdict, fontdict=fontdict, **kwargs)
    elif location == 14:
        t = ax.text(0.5, 1.0 + dist * aspect, label, horizontalalignment='center',
                verticalalignment='bottom', transform=ax.transAxes,
                bbox=boxdict, fontdict=fontdict, **kwargs)
    elif location == 15:
        t = ax.text(1 + dist, 0.5, label, horizontalalignment='left',
                verticalalignment='center_baseline', transform=ax.transAxes,
                bbox=boxdict, fontdict=fontdict, **kwargs)
    elif location == 16:
        t = ax.text(0.5, -dist * aspect, label, horizontalalignment='center',
                verticalalignment='top', transform=ax.transAxes,
                bbox=boxdict, fontdict=fontdict, **kwargs)
    elif location == 17:
        t = ax.text(- dist, 1.0 + dist * aspect, label,
                horizontalalignment='right', verticalalignment='bottom',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 18:
        t = ax.text(1.0 + dist, 1.0 + dist * aspect, label,
                horizontalalignment='left', verticalalignment='bottom',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 19:
        t = ax.text(1.0 + dist, 0.0 - dist * aspect, label,
                horizontalalignment='left', verticalalignment='top',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 20:
        t = ax.text(0.0 - dist, 0.0 - dist * aspect, label,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 21:
        t = ax.text(0.0 + dist, 0.5, label,
                horizontalalignment='left', verticalalignment='center_baseline',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 22:
        t = ax.text(0.5, 1.0 - dist * aspect, label,
                horizontalalignment='center', verticalalignment='top',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 23:
        t = ax.text(1.0 - dist, 0.5, label,
                horizontalalignment='right', verticalalignment='center_baseline',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    elif location == 24:
        t = ax.text(0.5, 0.0 + dist * aspect, label,
                horizontalalignment='center', verticalalignment='bottom',
                transform=ax.transAxes, bbox=boxdict,
                fontdict=fontdict, **kwargs)
    else:
        raise ValueError("Other corners not defined.")
    return t

def axes_from_axes(
        ax: Axes, n: int,
        extent: tp.Iterable = [0.2, 0.2, 0.6, 1.0],
        **kwargs) -> Axes:
    """Uses the location of an existing axes to create another axes in relative
    coordinates. IMPORTANT: Unlike ``inset_axes``, this function propagates
    ``*args`` and ``**kwargs`` to the ``pyplot.axes()`` function, which allows
    for the use of the projection ``keyword``.
    Parameters
    ----------
    ax : Axes
        Existing axes
    n : int
        label, necessary, because matplotlib will replace nonunique axes
    extent : list, optional
        new position in axes relative coordinates,
        by default [0.2, 0.2, 0.6, 1.0]
    Returns
    -------
    Axes
        New axes
    Notes
    -----
    DO NOT CHANGE THE INITIAL POSITION, this position works DO NOT CHANGE!
    :Author:
        Lucas Sawade (lsawade@princeton.edu)
    :Last Modified:
        2021.07.13 18.30
    """

    # Create new axes DO NOT CHANGE THIS INITIAL POSITION
    newax = plt.axes([0.0, 0.0, 0.25, 0.1], label=str(n), **kwargs)

    # Get new position
    ip = InsetPosition(ax, extent)

    # Set new position
    newax.set_axes_locator(ip)

    # return new axes
    return newax


def get_aspect(ax: Axes) -> float:
    """Returns the aspect ratio of an axes in a figure. This works around the
    problem of matplotlib's ``ax.get_aspect`` returning strings if set to
    'equal' for example
    Parameters
    ----------
    ax : Axes
        Matplotlib Axes object
    Returns
    -------
    float
        aspect ratio
    Notes
    -----
    :Author:
        Lucas Sawade (lsawade@princeton.edu)
    :Last Modified:
        2021.01.20 11.30
    """

    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()

    # Axis size on figure
    _, _, w, h = ax.get_position().bounds

    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)

    return disp_ratio


def set_default_color(COLOR):
    plt.rcParams['axes.edgecolor'] = COLOR
    plt.rcParams['text.color'] = COLOR
    plt.rcParams['axes.labelcolor'] = COLOR
    plt.rcParams['xtick.color'] = COLOR
    plt.rcParams['ytick.color'] = COLOR


def reset_mpl(gallery_conf, fname):
    """Function to set default look of the figures.
    """
    import matplotlib as mpl
    COLOR = 'k'
    mpl.rcParams["font.family"] = "monospace"
    mpl.rcParams["savefig.transparent"] = False
    mpl.rcParams["savefig.dpi"] = 300
    mpl.rcParams["savefig.format"] = 'svg'
    mpl.rcParams['axes.edgecolor'] = COLOR
    # mpl.rcParams['axes.backgroundcolor'] = 'w'
    mpl.rcParams['text.color'] = COLOR
    mpl.rcParams['axes.labelcolor'] = COLOR
    mpl.rcParams['xtick.color'] = COLOR
    mpl.rcParams['ytick.color'] = COLOR





def pick_colors_from_cmap(N: int, cmap: str = 'viridis') -> tp.List[tuple]:
    """Picks N uniformly distributed colors from a given colormap.

    Parameters
    ----------
    N : int
        Number of wanted colors
    cmap : str, optional
        name of the colormap to pick from, by default 'viridis'


    Returns
    -------
    List[tuple]
        List of color tuples.


    See Also
    --------
    lwsspy.plot.update_colorcycler.update_colorcycler : Updates the colors
        used in new lines/scatter points etc.

    """

    # Get cmap
    colormap = plt.get_cmap(cmap)

    # Pick
    colors = colormap(np.linspace(0, 1, N))

    return colors



import matplotlib.font_manager as fm
import matplotlib
import os
import glob

import matplotlib.ft2font as ft

FONTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fonts')


def updaterc(rebuild=False):
    """Updates the rcParams to something generic that looks ok good out of
    the box.

    Args:

        rebuild (bool):
            Rebuilds fontcache incase it needs it.

    Last modified: Lucas Sawade, 2020.09.15 01.00 (lsawade@princeton.edu)
    """

    add_fonts()

    params = {
        'font.family': 'sans-serif',
        'font.style':   'normal',
        'font.variant': 'normal',
        'font.weight':  'normal',
        'font.stretch': 'normal',
        'font.size':    12.0,
        'font.serif':     [
            'Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman',
            'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
            'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L',
            'Times', 'Palatino', 'Charter', 'serif'
        ],
        'font.sans-serif': [
            'Arial', 'Helvetica', 'DejaVu Sans', 'Bitstream Vera Sans',
            'Computer Modern Sans Serif', 'Lucida Grande', 'Verdana',
            'Geneva', 'Lucid', 'Avant Garde', 'sans-serif'
        ],
        'font.cursive':    [
            'Apple Chancery', 'Textile', 'Zapf Chancery', 'Sand', 'Script MT',
            'Felipa', 'Comic Neue', 'Comic Sans MS', 'cursive'
        ],
        'font.fantasy':    [
            'Chicago', 'Charcoal', 'Impact', 'Western', 'Humor Sans', 'xkcd',
            'fantasy'
        ],
        'font.monospace':  [
            'Roboto Mono', 'Monaco', 'DejaVu Sans Mono',
            'Bitstream Vera Sans Mono',  'Computer Modern Typewriter',
            'Andale Mono', 'Nimbus Mono L', 'Courier New', 'Courier', 'Fixed',
            'Terminal', 'monospace'
        ],
        'font.size': 12,
        # 'pdf.fonttype': 3,
        'figure.dpi': 140,
        'font.weight': 'normal',
        # 'pdf.fonttype': 42,
        # 'ps.fonttype': 42,
        # 'ps.useafm': True,
        # 'pdf.use14corefonts': True,
        'axes.unicode_minus': False,
        'axes.labelweight': 'normal',
        'axes.labelsize': 'small',
        'axes.titlesize': 'medium',
        'axes.linewidth': 1,
        'axes.grid': False,
        'grid.color': "k",
        'grid.linestyle': ":",
        'grid.alpha': 0.7,
        'xtick.labelsize': 'small',
        'xtick.direction': 'out',
        'xtick.top': True,  # draw label on the top
        'xtick.bottom': True,  # draw label on the bottom
        'xtick.minor.visible': True,
        'xtick.major.top': True,  # draw x axis top major ticks
        'xtick.major.bottom': True,  # draw x axis bottom major ticks
        'xtick.major.size': 4,  # draw x axis top major ticks
        'xtick.major.width': 1,  # draw x axis top major ticks
        'xtick.minor.top': True,  # draw x axis top minor ticks
        'xtick.minor.bottom': True,  # draw x axis bottom minor ticks
        'xtick.minor.width': 1,  # draw x axis top major ticks
        'xtick.minor.size': 2,  # draw x axis top major ticks
        'ytick.labelsize': 'small',
        'ytick.direction': 'out',
        'ytick.left': True,  # draw label on the top
        'ytick.right': True,  # draw label on the bottom
        'ytick.minor.visible': True,
        'ytick.major.left': True,  # draw x axis top major ticks
        'ytick.major.right': True,  # draw x axis bottom major ticks
        'ytick.major.size': 4,  # draw x axis top major ticks
        'ytick.major.width': 1,  # draw x axis top major ticks
        'ytick.minor.left': True,  # draw x axis top minor ticks
        'ytick.minor.right': True,  # draw x axis bottom minor ticks
        'ytick.minor.size': 2,  # draw x axis top major ticks
        'ytick.minor.width': 1,  # draw x axis top major ticks
        'legend.fancybox': False,
        'legend.frameon': True,
        'legend.loc': 'best',
        'legend.numpoints': 1,
        'legend.fontsize': 'small',
        'legend.framealpha': 1,
        'legend.scatterpoints': 3,
        'legend.edgecolor': 'inherit',
        'legend.facecolor': 'w',
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'sans',
        'mathtext.it': 'sans:italic',
        'mathtext.bf': 'sans:bold',
        'mathtext.cal': 'cursive',
        'mathtext.tt':  'monospace',
        'mathtext.default': 'it'
    }

    matplotlib.rcParams.update(params)


def add_fonts(verbose: bool = False):

    # Remove fontlist:
    for file in glob.glob('~/.matplotlib/font*.json'):
        os.remove(file)

    # Fonts
    fontfiles = glob.glob(os.path.join(FONTS, "*.tt?"))

    # for name, fname in fontdict.items():
    for fname in fontfiles:

        font = ft.FT2Font(fname)

        # Just to verify what kind of fonts are added verifiably
        if verbose:
            print(fname, "Scalable:", font.scalable)
            for style in ('Italic',
                          'Bold',
                          'Scalable',
                          'Fixed sizes',
                          'Fixed width',
                          'SFNT',
                          'Horizontal',
                          'Vertical',
                          'Kerning',
                          'Fast glyphs',
                          'Multiple masters',
                          'Glyph names',
                          'External stream'):
                bitpos = getattr(ft, style.replace(' ', '_').upper()) - 1
                print(f"{style+':':17}", bool(font.style_flags & (1 << bitpos)))

        # Actually adding the fonts
        fe = fm.ttfFontProperty(font)
        fm.fontManager.ttflist.insert(0, fe)

    # matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
    # def _add_Helvetica():

    #     # Check if Helvetica in system fonts
    #     from matplotlib import font_manager
    #     fonts = [os.path.basename(x).split(".")[0]
    #              for x in font_manager.findSystemFonts(
    #         fontpaths=None)]
    #     fonts.sort()
    #     # print(fonts)
    #     if "HelveticaNeue" in fonts:
    #         pass
    #     elif "Helvetica Neue" in fonts:
    #         pass
    #     elif "Helvetica" in fonts:
    #         return "Helvetica"
    #     else:
    #         font_file = os.path.join(
    #             os.path.dirname(__file__), 'fonts', 'HelveticaNeue.ttc')
    #         font_manager.fontManager.addfont(font_file)
    #     return "Helvetica Neue"

