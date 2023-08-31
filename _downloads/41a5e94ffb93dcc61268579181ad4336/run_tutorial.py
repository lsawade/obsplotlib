#!/bin/env python
"""
Full ``obsplotlib`` Tutorial
============================

The tutorial will go over the plotting functions in ``obsplotlib`` and how to
prepare your data to plot traces, station seismograms and full seismic sections
using ``obsplotlib``. Here and there the tutorial will digress into some
matplotlib details, to show how you could easily customize the plots to your
liking.


Loading all modules
-------------------

"""
# sphinx_gallery_thumbnail_number = 1
# sphinx_gallery_dummy_images = 1

import obspy
import obsplotlib.plot as opl
import matplotlib.pyplot as plt

# %%
#
# Loading data
# ------------
#
event = obspy.read_events("DATA/CMTSOLUTION")[0]
raw = obspy.read("DATA/observed/traces/*.sac")
inv = obspy.read_inventory("DATA/observed/station.xml")

# %%
# Before plotting anything let's get some information about the event and the
# station and process the data

# Get event latitude and longitude for geometry evaluation
event_time = event.preferred_origin().time
event_latitude = event.preferred_origin().latitude
event_longitude = event.preferred_origin().longitude
event_depth = event.preferred_origin().depth  # in meters
event_name = 'C' + event.preferred_origin().resource_id.id.split('/')[-2]

# %%
# Attach the event station geometry to the traces, important for rotation to RTZ

opl.attach_geometry(raw, event_latitude=event_latitude,
                    event_longitude=event_longitude, inv=inv)

# %%
# Processing the data very generically
bandpass = [30, 200]
obs = opl.process(raw, inv=inv, remove_response=True, bandpass=bandpass)

# %%
# Trace
# -----
#
# First let's inspect the trace plotting function. It is in a way a wrapper
# around the matplotlib plot function but with some added functionality to grab
# info from the stats object and plot it in a "nice" way.

# Select a trace
network_str, station_str, component_str = "II", "BFO", "Z"
tr = obs.select(network=network_str, station=station_str,
                component=component_str)[0]

plt.figure()
ax = opl.trace(tr, plot_labels=True, lw=0.5)
plt.show(block=False)

# %%
# Since we have both station and event information in the stats object we can
# add a header to the figure to be a little more explicit.

plt.figure(figsize=(8, 3))
ax = opl.trace(tr, plot_labels=False, origin_time=event_time, lw=0.5)

header_dict = dict(
    station=f'{tr.id}',
    station_latitude=tr.stats.latitude,
    station_longitude=tr.stats.longitude,
    station_azimuth=tr.stats.azimuth,
    station_back_azimuth=tr.stats.back_azimuth,
    station_distance_in_degree=tr.stats.distance,
    event=event_name,
    event_latitude=event_latitude,
    event_longitude=event_longitude,
    event_depth_in_km=event_depth/1000.0,
    event_time=event_time,
    add_newline_station=True,
    add_newline_event=True,
    bandpass=bandpass,
    fontsize='medium'
)

opl.add_header(ax, **header_dict)
plt.subplots_adjust(left=0.05, right=0.95, top=0.725, bottom=0.15)
plt.show(block=False)

# %%
# In the example above we are providing all possible arguments to the
# ``add_header`` function just for show. Depending on whether they are provided
# they will be added to the header or not. The header is a simple text object
# and all font related arguments are passed through ``plot_label`` and to
# ``plt.text()``. The ``add_newline_station`` and ``add_newline_event`` arguments
# Simply add a newline and a space after the station name and event name.
#
# .. note::
#
#     Digression: At this point you probably already noticed how I'm using
#     a monospace font. You may adjust this to your liking by changing the
#     ``plt.rcParams["font.family"]`` parameter, e.g
#
#     .. code:: python
#
#         plt.rcParams["font.family"] = "Arial"
#
#     Monospace is a personal preference of mine, because it makes it easier to
#     align the header and the labels. But it is not the most beautiful font.
#     Especially, if you are comparing traces and plot labels that contain
#     numbers, it is simpler to compare numbers if they are aligned. Anywho
#     I will enable it for the next section before switching back to monospace.
#
# Station
# -------
#
# The next function is the station function. Instead of plotting a single trace
# it will plot a set of components in a single figure. It's a wrapper around the
# plot trace function, so most arguments are parsed to the trace function. The
# The components are defined by a keyword argument. So you may use ``ZRT``
# ``NEZ`` or ``123`` or just two, ``RT``, for example.

# Switching font to Arial
plt.rcParams["font.family"] = "Arial"

# Get station from observed trace
st = obs.select(network="II", station="BFO")

# Plot the station
plt.figure(figsize=(8, 5))
axes = opl.station(st, components='ZRT', lw=0.5)

# If dissatisfied with legend fontsize and position? Just recreate it using
# the first axes object.
axes[0].legend(frameon=False, loc='lower right', ncol=3, fontsize='small',
               bbox_to_anchor=(1.0, 1.0))

# Add the header with a bit more distance to make room for the legend outside
# the  axes
opl.add_header(axes[0], **header_dict, dist=0.075)

# Slightly adjust the plots to make the fit nicely into the figure
plt.subplots_adjust(left=0.075, right=0.925, top=0.775, bottom=0.15)
plt.show(block=False)

# %%
# Putting the legend outside the axes is nice when we are plotting multiple
# traces to compare them. But if we are only plotting a single stream, it is
# nicer to plot the legend inside the axes. Because it removes some unused
# white space.
#
# Section
# -------
#
# Plotting a section should be simple. And obspy does make it fairly easy, but
# the moment you want to plot a section with multiple components, or align
# traces it becomes fairly complicated. ``obsplotlib`` is trying to streamline
# these processes, by using some function to add properties to the stats object
# of the traces and then using these properties to plot the section. We actually
# already did this earlier in the tutorial when we used
# ``opl.attach_geometry(...)`` to attach the station coordinates to the stats.
# If you don't have a stationxml file, you can use opl.attach_geometry(...)
# after attaching station coordinates to the traces' stats objects, or manually
# add distance, (and optionally, azimuth, and back_azimuth for labels)
# to the stats object. To save space in the section, traces are not plotting by
# their actual distance, but one by one with a label that has the distance.

# Switching back to monospace
plt.rcParams["font.family"] = "monospace"

plt.figure(figsize=(8, 10))
opl.section(obs, lw=0.5, comp='T')
plt.legend(frameon=False, loc='upper right', ncol=3, fontsize='small')
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)
plt.show(block=False)


# %%
# Plotting the same section but with axis limits from 300 seconds to 1500 seconds
# after the event
starttime = event_time + 600
endtime = event_time + 3600
limits = [starttime, endtime]

plt.figure(figsize=(8, 10))
opl.section(obs, limits=limits, lw=0.5)
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)
plt.show(block=False)


# %%
# Plotting the same section but with the origin time defined so that we get
# Time since event origin on the x-axis.

plt.figure(figsize=(8, 10))
opl.section(obs, origin_time=event_time, lw=0.5)
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)
plt.show(block=False)

# %%
# Again plotting the same section but now with the origin time defined so that we
# get Time since event origin on the x-axis and with axis limits from 300 seconds
# to 1500 seconds after the event

limits = [600, 3600]
plt.figure(figsize=(8, 10))
opl.section(obs, origin_time=event_time, limits=limits, lw=0.5)
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)
plt.show(block=False)

# %%
# Next we are going to plot an aligned sections. To do this each trace must have
# a obspy.Trace.stats.traveltime parameter. This can be done using the
# add_traveltime function or manually using your own function.

obs_filtered = opl.add_traveltime(obs, phase='love', orbit=1,
                                  return_filtered=True, vlove=6.5)

# %%
# Note that the add traveltime function uses the TauPy model by default for
# body waves and fixed velocity for surface waves. The traveltimes are then
# computed using the distance parameter in the stats object.

# %%
plt.figure(figsize=(8, 10))
opl.section(obs_filtered, origin_time=event_time, lw=0.5, align=True,
            comp='T')
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)
plt.show(block=False)


# %%
# Now this does not make a lot of sense since the traces are not aligned at the
# start. Conveniently we can set the limits parameter to only plot a certain
# time range before and after the arrival times.

limits = [-500, 500]

plt.figure(figsize=(6, 10))
opl.section(obs_filtered, origin_time=event_time, limits=limits, lw=0.5,
            align=True, comp='T')
plt.subplots_adjust(left=0.25, right=0.8, top=0.95, bottom=0.05)
plt.show(block=False)

# %%
# Let's also plot a section aligned to the P arrival times.

obs_filtered = opl.add_traveltime(obs, event_depth_in_m=event_depth, phase='P',
                                  origin_time=event_time,
                                  return_filtered=True, vlove=3.7)
limits = [-100, 150]

plt.figure(figsize=(8, 6))
ax, ax2 = opl.section(obs_filtered, origin_time=event_time, limits=limits, lw=0.5,
                      align=True)
ax.plot([0, 0], [-0.5, len(obs_filtered) + 0.5], 'k--', lw=0.5)
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.1)
plt.show(block=False)

# %%
# Note that we have fewer traces here because some land in the Pwave shadow zone
# and are not recorded but seismographs. Also note, that ax, and ax2 give you
# axes to the left and right yaxes. ax and ax2 ticks actually set the left and
# right y axes labels. So far we have only plotted a single
# component (Z) in the section. ``obsplotlib`` also has a function to plot
# multiple components in a single section. This is done using the
# ``opl.section_multiple_comp`` function. This function takes the same arguments
# But insted of being a single letter string, the components argument is a
# string with all components to be plotted in order.

plt.figure(figsize=(9, 5))
axes = opl.section_multiple_comp(obs_filtered, origin_time=event_time,
                                 limits=limits, lw=0.5, align=True,
                                 components="ZRT")
for ax, _ in axes:
    ax.plot([0, 0], [-0.5, len(obs_filtered) + 0.5], 'k--', lw=0.5)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05, wspace=0.75)
plt.show(block=False)

# %%
# One main difference is that the section multiple components will find an
# absmax to normalize across all streams and traces. This can be overwritten by
# absmax parameter which can be manually set.
#
#
# Trace comparison
# ----------------
#
# So far we have only really looked at a single set of traces. Very often in
# seismology however we want to look at trace comparisons. And sometimes
# directly looks at measurements on traces, or windows. Let's load a second
# set of traces to compare our observed data too.

# Read traces and station info
synraw = obspy.read("DATA/synthetic/traces/*.sac")
syninv = obspy.read_inventory("DATA/synthetic/station.xml")

# %%
# Just like with the observed data we are attach geometry for rotation
opl.attach_geometry(synraw, event_latitude=event_latitude,
                    event_longitude=event_longitude, inv=syninv)


# %%
# Since we want to process both synthetics and observed the same fashion,
# We have to resample the traces in addition to the basic processing.

starttime = event_time
npts = 10800
sampling_rate_in_hz = 1
bandpass = [50, 500]
obs = opl.process(raw, inv=inv, remove_response=True, bandpass=bandpass,
                  starttime=starttime, npts=npts,
                  sampling_rate_in_hz=1)

syn = opl.process(synraw, inv=inv, remove_response=False, bandpass=bandpass,
                  starttime=starttime, npts=npts,
                  sampling_rate_in_hz=1)

# %%
# Once both are processed we can plot them with

obstr = obs.select(network=network_str, station=station_str,
                   component=component_str)[0]
syntr = syn.select(network=network_str, station=station_str,
                   component=component_str)[0]

plt.figure()
ax = opl.trace([obstr, syntr], labels=['Observed', 'GLAD-M25'],
               origin_time=event_time, lw=0.75)

# Just reusing the header dict from earlier
header_dict['station'] = obstr.id
header_dict['bandpass'] = bandpass

opl.add_header(ax, **header_dict)
plt.subplots_adjust(left=0.05, right=0.95, top=0.8, bottom=0.125)
plt.show(block=False)

# %%
# Repeat to plot a station

# Get station from observed trace
obs_st = obs.select(network="II", station="BFO")
syn_st = syn.select(network="II", station="BFO")

# Plot the station
plt.figure(figsize=(8, 5))
axes = opl.station([obs_st, syn_st], components='ZRT', lw=0.5,
                   labels=['Observed', 'GLAD-M25'], nooffset=True)

# If dissatisfied with legend fontsize and position? Just recreate it using
# the first axes object.
axes[0].legend(frameon=False, loc='lower right', ncol=3, fontsize='small',
               bbox_to_anchor=(1.0, 1.0))

# Add the header with a bit more distance to make room for the legend outside
# the  axes
opl.add_header(axes[0], **header_dict, dist=0.075)

# Slightly adjust the plots to make the fit nicely into the figure
plt.subplots_adjust(left=0.075, right=0.925, top=0.775, bottom=0.15)
plt.show(block=False)

# %%
# For the section, we need do a couple more things. The set of traces in these
# do not perfectly overlapping.

streams = opl.select_intersection([obs, syn], components='ZRT')

# %%
# Now that we have selected only the traces that are in all streams, we can
# plot a section with the traces.

plt.figure(figsize=(8, 10))
opl.section(streams, origin_time=event_time, lw=0.5,
            labels=['Observed', 'GLAD-M25'])
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)
plt.show(block=False)
