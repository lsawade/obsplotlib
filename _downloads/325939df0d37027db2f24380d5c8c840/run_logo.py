#!/bin/env python
"""
Plotting the logo and the favicon
=================================

This example will go over the creation of the logo from the data in the example
directory.


Loading all modules
-------------------

"""
# sphinx_gallery_dummy_images = 1

import numpy as np
import obspy
import obsplotlib.plot as opl
import matplotlib.pyplot as plt

# %%
#
# Loading data
# ------------
#
event = obspy.read_events("DATA/CMTSOLUTION")[0]
raw = obspy.read("DATA/observed/traces/II.BFO.*.sac")
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
# Processing data
# ---------------
#
# Here we're just taking the data from the tutorial and are playing with the
# taper.
#
bandpass = [10, 50]
obs = opl.process(raw, inv=inv, remove_response=True, bandpass=bandpass,
                  starttime=event_time+800, npts=750, sampling_rate_in_hz=5
                  )

# Select a trace
network_str, station_str, component_str = "II", "BFO", "Z"
tr = obs.select(network=network_str, station=station_str,
                component=component_str)[0]

# %%
# Finally we perform some fixes to really isolate the wavelet

# Tapering multiple time to attenuate boundaries
for i in range(10):
    tr.taper(type='cosine', max_percentage=0.4)

# Adding zero padding
tr.data = np.hstack((np.zeros(150), tr.data[:-75]))

# Offset the trace to be at the height of the 'T'
tr.data += 0.0000015


# %%
# Plot Logo
# ---------
#
# We want the logo to be visible on both dark and light backgrounds. As it turns
# out, soft orange and soft blue are ideal colors for that

orange = np.array([232, 152, 59])/255
blue = np.array([65, 140, 216])/255

plt.figure(figsize=(2.0, 0.5))
ax = opl.trace(tr, plot_labels=False, lw=1.4, origin_time=event_time,
               solid_joinstyle='round', solid_capstyle='round', colors=[orange,],
               zorder=10)

# Remove bottom axis
ax.spines['bottom'].set_visible(False)

# Remove ticks
ax.tick_params(axis='x', labelbottom=False, bottom=False)

# Remove legend
ax.get_legend().set_visible(False)

# Remove xlabel
plt.xlabel('')

# Add label for
opl.plot_label(ax, 'Obsplotlib', location=13,
               dist=-0.41, box=False, fontsize='x-large', color=blue)
plt.subplots_adjust(left=0.35, top=0.95, bottom=0.05, right=1.0)
plt.show(block=False)

# To store as SVG
plt.savefig('logo.svg', transparent=True)


# %%
# Plot Favicon
# ------------
#
# We want the logo to be visible on both dark and light backgrounds. As it turns
# out, soft orange and soft blue are ideal colors for that

orange = np.array([232, 152, 59])/255
blue = np.array([65, 140, 216])/255

plt.figure(figsize=(1.0, 1.0))


ax = opl.trace(tr, plot_labels=False, lw=4.0, origin_time=event_time,
               solid_joinstyle='round', solid_capstyle='round', colors=[orange,],
               zorder=10)
ax = opl.trace(tr, plot_labels=False, lw=4.0, origin_time=event_time-2,
               solid_joinstyle='round', solid_capstyle='round', colors=[blue,],
               zorder=10)

# Remove bottom axis
ax.spines['bottom'].set_visible(False)

# Remove ticks
ax.tick_params(axis='x', labelbottom=False, bottom=False)

# Remove legend
ax.get_legend().set_visible(False)

# Remove xlabel
plt.xlabel('')
plt.xlim(867, 905)
plt.subplots_adjust(left=0.0, top=1.0, bottom=0.0, right=1.0)
plt.show(block=False)

# To store as SVG
# plt.savefig('favicon.ico', format='PNG', transparent=True, dpi=32)
