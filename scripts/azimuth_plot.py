# %%
import numpy as np
import matplotlib.pyplot as plt


# %%
# Create a function that generates a set of 20 random angles between 0, 2pi
# but make sure that they are somwhat evenly distributed
def generate_angles(N):
    import numpy as np

    # Generate 20 random angles
    angles = np.random.rand(N) * 2 * np.pi

    # Sort the angles
    angles = np.sort(angles)

    return angles


N = 100
angles = generate_angles(N)


# %%


# Create a function that generates a set of N timeseries, each with a random
# choice of either a gaussian, a trapezoid or a triangular shape. The signals
# should not exceed 10% of the time series since the start, and the time series
# should have length NT
def generate_timeseries(N: int, NT: int, factor=1.0):
    import numpy as np

    # Generate a random set of shapes
    shapes = np.random.choice(["gaussian"], N)

    # Generate a random set of centroid times
    t0 = np.random.rand(N) * (NT * 0.05) + 0.1 * NT

    # Generate a random set of M0 values
    A = (np.random.rand(N) * 0.25 + 1) * factor

    # Generate a random set of hdur values
    alpha = 1.628
    sigma = (np.random.rand(N) * 5 + 15) / alpha

    # Generate a time vector
    t = np.arange(NT)

    # Generate an empty array to store the signals
    signals = np.zeros((N, NT))

    # Loop over the number of signals
    for i in range(N):

        # Check the shape
        if shapes[i] == "gaussian":

            # Exponent for the Gaussian
            exponent = -(((t - t0[i]) / sigma[i]) ** 2)

            # Are under the Gaussen -> M0
            gaussian = A[i] / (np.sqrt(np.pi) * sigma[i]) * np.exp(exponent)

            # Store the signal
            signals[i] = gaussian

        elif shapes[i] == "trapezoid":

            # Startpoint
            startpoint = t0[i] - sigma[i]

            # Midpoint
            midpoint = t0[i]

            # Endpoint
            endpoint = t0[i] + sigma[i]

            # total area under triangle has to be 1
            trapezoid = np.zeros(NT)
            trapezoid[(t >= startpoint) & (t <= midpoint)] = A[i] / sigma[i]
            trapezoid[(t > midpoint) & (t <= endpoint)] = (
                A[i]
                / sigma[i]
                * (1 - (t[(t > midpoint) & (t <= endpoint)] - midpoint) / sigma[i])
            )

            # Store the signal
            signals[i] = trapezoid

        elif shapes[i] == "triangular":

            # Startpoint
            startpoint = t0[i] - sigma[i]

            # Midpoint
            midpoint = t0[i]

            # Endpoint
            endpoint = t0[i] + sigma[i]

            # total area under triangle has to be 1
            triangle = np.zeros(NT)

            # First half
            triangle[(t >= startpoint) & (t <= midpoint)] = A[i] * (
                t[(t >= startpoint) & (t <= midpoint)] - startpoint
            )
            # Second half
            triangle[(t > midpoint) & (t <= endpoint)] = A[i] * (
                1 - (midpointt[(t > midpoint) & (t <= endpoint)])
            )

            # Store the signal
            signals[i] = triangle

    return signals


signals = generate_timeseries(N, 1000) + generate_timeseries(N, 1000, factor=-1.0)

# %%


# We generated the angles. Now we want to bin the angles into 5 degree bins
# and select one angle from each bin but ignore bins that do not contain angles
def bin_angles(angles, binwidth: float):
    import numpy as np

    # Create a list to store the new angles
    new_angles = []
    new_indeces = []

    # Loop over the bins
    for i in np.arange(0, 360 + binwidth, binwidth):

        # Get the angles in the bin
        bin_angles = angles[
            (angles >= np.deg2rad(i)) & (angles < np.deg2rad(i + binwidth))
        ]

        if len(bin_angles) == 0:
            continue

        # for later I also need the indices
        indices = np.where(
            (angles >= np.deg2rad(i)) & (angles < np.deg2rad(i + binwidth))
        )

        # append
        new_angles.append(np.mean(bin_angles))
        new_indeces.append(indices[0][0])
    return new_angles, new_indeces


binwidth = 10
new_angles, new_indeces = bin_angles(angles, 10)
new_signals = signals[new_indeces]

# %%

# Now I want to create a figure that creates one axes for each signal, and
# and plots the signal in the axes. But the axes should be arranged with respect to azimuth

import obsplotlib.plot as opl


plt.figure(figsize=(10, 10))
mainax = plt.gca()
mainax.axis("off")

bins = np.arange(0, 360 + binwidth, binwidth)

for i in range(len(new_angles)):

    # length of the axes
    width = 0.3
    height = 0.05

    # Distance of the axes from the center
    r = 0.3
    stretch_height = 1.0
    relative_stretch = 1.0

    # Angle of the axes
    az = new_angles[i]

    # Use azimuth to get x,y. Here azimuth is with respect to North and
    # clockkwise
    x = r * np.sin(az) + 0.5
    y = stretch_height * (1 + relative_stretch * np.cos(az) ** 2) * r * np.cos(az) + 0.5

    # plot bin edges
    for _i in bins[:-1]:
        mainax.plot(
            [0.5, 0.5 + r * np.sin(np.deg2rad(_i))],
            [
                0.5,
                0.5
                + stretch_height
                * (1 + relative_stretch * np.cos(np.deg2rad(_i)) ** 2)
                * r
                * np.cos(np.deg2rad(_i)),
            ],
            c="lightgray",
            lw=0.1,
            ls="-",
            clip_on=False,
        )

    print(x, y, az)

    # If the azimuth is larger than pi the axis is on the left side
    axis_left = az >= np.pi
    ylabel_location = "right" if axis_left else "left"
    yrotation = 0 if axis_left else 0

    # ADjust the left axes to match the width
    if axis_left:
        x = x - width - width * 0.1
    else:
        x = x + width * 0.1

    # Create extent
    extent = [x, y - height / 2, width, height]
    print(extent)

    # Create axes
    ax = opl.axes_from_axes(
        mainax,
        9809834 + i,
        extent=extent,
    )
    # Remove all ticks and labels
    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labeltop=False,
        labelleft=False,
        labelright=False,
    )
    # Remove spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    # ax.spines['left'].set_visible(False)

    # Remove axis background
    ax.patch.set_visible(False)

    # Set the ylabel
    ax.yaxis.set_label_position(ylabel_location)
    ax.plot(new_signals[i, 50:250])
    ax.set_ylabel(
        f"{az*180/np.pi:.0f}",
        rotation=yrotation,
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=10,
        fontsize="small",
    )
    mainax.scatter(
        x,
        y,
        s=5,
        c=az,
        marker="o",
        cmap="viridis",
        vmin=0,
        vmax=2 * np.pi,
        clip_on=False,
        zorder=20,
    )
    mainax.set_xlim(0, 1)
    mainax.set_ylim(0, 1)
