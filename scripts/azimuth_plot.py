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
def generate_timeseries(angles, NT: int, factor=1.0, strike=0.0):
    import numpy as np

    # Fault parameters
    L = 10  # length of fault
    vR = 2.0  # rupture velocity
    v = 1.2 * vR

    # Number of signals
    N = len(angles)

    # Generate a random set of shapes
    shapes = np.random.choice(["gaussian"], N)

    # Generate a random set of centroid times
    t0 = np.random.rand(N) * (NT * 0.05) + 0.1 * NT

    # Generate a random set of M0 values
    A = (np.random.rand(N) * 0.25 + 1) * factor

    # Generate a random set of hdur values
    alpha = 1.628

    # Strike dependent hdur
    hdur = L * (1 / vR - np.cos((angles - np.deg2rad(strike))) / v)

    sigma = (hdur) / alpha

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
                1 - (midpoint[(t > midpoint) & (t <= endpoint)])
            )

            # Store the signal
            signals[i] = triangle

    return signals


signals = generate_timeseries(
    angles, 1000, strike=45.0
)  # + generate_timeseries(N, 1000, factor=-1.0)

# %%


# We generated the angles. Now we want to bin the angles into 5 degree bins
# and select one angle from each bin but ignore bins that do not contain angles
def bin_angles(angles, dy: float = 0.1):
    import numpy as np

    y = np.clip(np.arange(1, -1 - dy, -dy), -1, 1)

    theta = np.arccos(y)
    x = np.sin(theta)

    # Make sure that we wave bins for both sides, positive and negative x
    x = np.concatenate([x[:-1], -x[::-1]])
    y = np.concatenate([y[:-1], y[::-1]])
    theta = np.concatenate([theta[:-1], np.pi + theta])

    # Initialize the new angles, and indeces
    new_angles = []
    new_indeces = []

    # Loop over the bins
    for i in range(len(theta) - 1):

        # Get the angles in the bin
        bin_angles = angles[(angles >= theta[i]) & (angles < theta[i + 1])]

        if len(bin_angles) == 0:
            continue

        # for later I also need the indices
        indices = np.where((angles >= theta[i]) & (angles < theta[i + 1]))

        # Choose the angle that is the closest to the center of the bin
        center = (theta[i] + theta[i + 1]) / 2
        index = np.argmin(np.abs(bin_angles - center))

        # append
        new_angles.append(bin_angles[index])
        new_indeces.append(indices[0][index])

    return (theta, x, y), (new_angles, new_indeces)


dy = 0.1
(theta_bin, x_bin, y_bin), (new_angles, new_indeces) = bin_angles(angles, dy=dy)
new_signals = signals[new_indeces]

# %%

# Now I want to create a figure that creates one axes for each signal, and
# and plots the signal in the axes. But the axes should be arranged with respect to azimuth

import obsplotlib.plot as opl


plt.figure(figsize=(10, 10))
mainax = plt.gca()
mainax.axis("off")

for i in range(len(new_angles)):

    # Height ration

    # Height of the axes as a function of stretch
    stretch_height = 2.0
    r = 0.5 / stretch_height
    height = 0.5 * dy

    # Distance of the axes from the center
    # r = 0.2
    total_width = 0.5 - r
    percentage_offset = 0.1
    width = total_width * (1 - percentage_offset)
    width_offset = total_width * percentage_offset

    # Height of the axes as a function of stretch
    # stretch_height = 2.0
    # height = stretch_height * dy * r

    # Angle of the axes
    az = new_angles[i]

    # Use azimuth to get x,y. Here azimuth is with respect to North and
    # clockkwise
    x = r * np.sin(az) + 0.5
    y = stretch_height * r * np.cos(az) + 0.5

    # plot bin edges
    for _i in theta_bin[:-1]:
        mainax.plot(
            [0.5, 0.5 + r * np.sin(_i)],
            [0.5, 0.5 + stretch_height * r * np.cos(_i)],
            c="lightgray",
            lw=0.1,
            ls="-",
            clip_on=False,
        )

    # If the azimuth is larger than pi the axis is on the left side
    axis_left = az >= np.pi
    ylabel_location = "right" if axis_left else "left"
    yrotation = 0 if axis_left else 0

    # ADjust the left axes to match the width
    if axis_left:
        x = x - width - width_offset
    else:
        x = x + width_offset

    # Create extent
    extent = [x, y - height / 2, width, height]

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

# mainax.set_aspect("equal")
plt.subplots_adjust(left=0.0, right=1.0, top=0.95, bottom=0.05)


# %%
plt.ion()
plt.figure()
plt.plot(x, y, "o")
plt.gca().set_aspect("equal")
plt.show()
