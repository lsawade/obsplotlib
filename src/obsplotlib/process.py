# Some functions
import obspy


def process(st, inv: obspy.Inventory | None = None,
            remove_response=False, to_RTZ=True,
            starttime=None, npts=None,
            sampling_rate_in_hz: float = 1,
            bandpass=[200, 500]):

    out = st.copy()

    # Basics
    out.detrend('demean')
    out.detrend('linear')
    out.taper(max_percentage=0.05, type='cosine')

    # Remove response
    if inv is not None and remove_response:
        out.remove_response(inventory=inv, output="DISP",
                            pre_filt=(0.001, 0.005, 0.1, 0.2))

    # Rotate to NEZ (for standardization)
    if inv is not None:
        out.rotate('->ZNE', inventory=inv)

    # Bandpass filter traces given the bandpass parameters
    out.filter('bandpass',
               freqmin=1 / bandpass[1], freqmax=1 / bandpass[0],
               corners=2, zerophase=True)

    # Rotate to RTZ
    if to_RTZ:
        out.rotate('NE->RT', inventory=inv)

    # Interpolate to standardizeif
    if isinstance(starttime, obspy.UTCDateTime) \
            and isinstance(npts, int):
        out.interpolate(starttime=starttime, npts=npts,
                        sampling_rate=sampling_rate_in_hz)

    out.taper(max_percentage=0.05, type='cosine')

    return out
