#!/bin/env python
"""
Aligned Station Section
=======================

The tutorial will go over the plotting an aligned section of waveforms
using built-in plotting tools.

Loading all modules
-------------------


"""
# sphinx_gallery_thumbnail_number = 1
# sphinx_gallery_dummy_images = 1
# %%
from obsplotlib.seismogram import station
import obspy
import obsplotlib.plot as opl
import matplotlib.pyplot as plt

event = obspy.read_events("DATA/CMTSOLUTION")[0]
raw = obspy.read("DATA/observed/traces/*.sac")
inv = obspy.read_inventory("DATA/observed/station.xml")

# %%
# Get event latitude and longitude for geometry evaluation
event_latitude = event.preferred_origin().latitude
event_longitude = event.preferred_origin().longitude

# %%
# Attach the event station geometry to the traces, important for rotation to RTZ
opl.attach_geometry(raw, event_latitude=event_latitude,
                    event_longitude=event_longitude, inv=inv)

# %%
# Processing the data very generically
obs = opl.process(raw, inv=inv, remove_response=True, bandpass=[20, 200])

# %%
# Plotting a section should be simple. Here the obs object should already contain
# distance parameters for each trace. If not, use opl.attach_geometry(...)
# If you don't have a stationxml file, you can use opl.attach_geometry(...)
# after attaching station coordinates to the traces' stats objects.

plt.figure()
opl.section(obs)
plt.show(block=False)

# %%
# Plot single trace
tr = obs.select(network="II", station="BFO", component='Z')[0]
plt.figure()
ax = opl.trace(tr, plot_labels=False)
plt.show(block=False)


# %%
# Plot stations
# %%
st = obs.select(network="II", station="BFO")
plt.figure()
ax = opl.station(st)
plt.show(block=False)
