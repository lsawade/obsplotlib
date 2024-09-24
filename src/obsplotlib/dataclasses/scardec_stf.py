from __future__ import annotations
import os
from glob import glob
import obspy
import numpy as np
from dataclasses import dataclass
import obspy.core.event.event


@dataclass
class SCARDECSTF:
    """Inside each earthquake directory, two files are provided, for the average STF (file fctmoysource_YYYYMMDD_HHMMSS_Name) and for the optimal STF (file fctoptsource_YYYYMMDD_HHMMSS_Name)

     These two STF files have the same format:

    1st line: YYYY MM DD HH MM SS'.0' Latitude Longitude [origin time and epicentral location from NEIC]
    2nd line: Depth(km) M0(N.m) Mw strike1(°) dip1(°) rake1(°) strike2(°) dip2(°) rake2(°) [all from SCARDEC]
    All the other lines are the temporal STF, with format: time(s), moment rate(N.m/s)
    """

    origin: obspy.UTCDateTime
    latitude: float
    longitude: float
    depth_in_km: float
    M0: float
    Mw: float
    strike1: float
    dip1: float
    rake1: float
    strike2: float
    dip2: float
    rake2: float
    time: np.ndarray
    moment_rate: np.ndarray
    region: str

    @classmethod
    def fromfile(cls, filename):

        # Get region from filename
        region = " ".join(os.path.basename(filename).split("_")[3:])

        with open(filename, "r") as filename:

            lines = filename.readlines()

        line1 = lines[0].split()
        line2 = lines[1].split()

        origin = obspy.UTCDateTime(
            int(line1[0]),
            int(line1[1]),
            int(line1[2]),
            int(line1[3]),
            int(line1[4]),
            float(line1[5]),
        )
        latitude = float(line1[6])
        longitude = float(line1[7])

        depth_in_km = float(line2[0])
        M0 = float(line2[1])
        Mw = float(line2[2])
        strike1 = int(line2[3])
        dip1 = int(line2[4])
        rake1 = int(line2[5])
        strike2 = int(line2[6])
        dip2 = int(line2[7])
        rake2 = int(line2[8])

        # Now get STF
        time = []
        moment_rate = []
        for line in lines[2:]:
            t, m = line.split()
            time.append(float(t))
            moment_rate.append(float(m))

        # Convert to numpy arrays
        time = np.array(time)
        moment_rate = np.array(moment_rate)

        # Create the object
        return cls(
            origin,
            latitude,
            longitude,
            depth_in_km,
            M0,
            Mw,
            strike1,
            dip1,
            rake1,
            strike2,
            dip2,
            rake2,
            time,
            moment_rate,
            region,
        )

    @classmethod
    def fromdir(cls, dirname, stftype="optimal"):
        if stftype == "optimal":
            return cls.fromfile(glob(os.path.join(dirname, "fctoptsource*"))[0])
        elif stftype == "average":
            return cls.fromfile(glob(os.path.join(dirname, "fctmoysource*"))[0])
        else:
            raise ValueError("stftype must be 'optimal' or 'average'")
