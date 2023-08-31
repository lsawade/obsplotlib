import obspy
import numpy as np
import typing as tp


def check_traces(tr1: obspy.Trace, tr2: obspy.Trace):
    if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
        raise ValueError('Sampling rates must be the same')

    if tr1.stats.starttime != tr2.stats.starttime:
        raise ValueError('starttime must be the same')

    if tr1.stats.endtime != tr2.stats.endtime:
        raise ValueError('endtime must be the same')

    return True


def L2(tr1: obspy.Trace, tr2: obspy.Trace, normalize: bool = True) -> float:
    """Returns the L2 norm between two Traces."""

    # Raises error if traces are not compatible
    check_traces(tr1, tr2)

    # Return the L2 Norm
    misfit = np.sum((tr1.data-tr2.data)**2)

    # How to normalize the traces. Either number of samples or by first trace
    if normalize:
        norm = 1.0 / np.sum(tr1.data**2)
    else:
        norm = 0.5 * tr1.stats.delta / len(tr1.data)

    return float(norm * misfit)


def X(tr1: obspy.Trace, tr2: obspy.Trace, normalize: bool = True) \
        -> tp.Tuple[float, int]:
    return xcorr(tr1.data, tr2.data)


def xcorr(d, s) -> tp.Tuple[float, int]:
    """Compute cross-correlation between two signals.

    Parameters
    ----------
    d : ndarray
        signal to correlate
    s : ndarray
        signal to correlate with

    Returns
    -------
    tp.Tuple[float, int]
        max_cc_value, time_shift.
        The time shift is the positive if d is delayed. And the timeshift is
        negative if s is delayed. In other words if d and s are data and
        synthetics, your synthetics are late if the timeshit is negative.

    """
    cc = np.correlate(d, s, mode="full")
    time_shift = cc.argmax() - len(d) + 1
    # Normalized cross correlation.
    max_cc_value = cc.max() / np.sqrt((s ** 2).sum() * (d ** 2).sum())

    return max_cc_value, time_shift


def Xratio(tr1: obspy.Trace, tr2: obspy.Trace, normalize: bool = True):
    """Correlation ratio between two signals. -> (trace1*trace2)/(trace2**2).
    This provides a measure of how well the two signals correlate.

    Parameters
    ----------
    tr1 : obspy.Trace
        signal to correlate
    tr2 : obspy.Trace
        signal to correlate with
    normalize : bool, optional
        Normalize the cross-correlation, by default True

    Useful for comparing the amplitude of two signals after cross correlation
    time correction. The ratio is 1 if the two signals are identical. If
    trace 1 is larger than trace 2, the ratio is > 1. If trace 1 is smaller
    the ratio is < 1. Can be used for amplitude correction. See some Ekström
    paper related to surfacee wave measurements and station corrections.

    Returns
    -------
    float
        The ratio of the cross-correlation to the maximum cross-correlation
        value.

    """
    return np.sum(tr1.data * tr2.data)/np.sum(tr2.data ** 2)


def correct_window_index(istart, iend, nshift, npts):
    """Correct the window index based on cross-correlation shift

    Parameters
    ----------
    istart : int
        start index
    iend : int
        end index
    nshift : int
        shift in N samples
    npts : int
        Length of window

    Returns
    -------
    Tuple
        indeces

    Raises
    ------
    ValueError
        If resulting windows arent the same length? I don't get this

    """
    istart_d = max(1, istart + nshift)
    iend_d = min(npts, iend + nshift)
    istart_s = max(1, istart_d - nshift)
    iend_s = min(npts, iend_d - nshift)
    if (iend_d - istart_d) != (iend_s - istart_s):
        raise ValueError("After correction, window length not the same: "
                         "[%d, %d] and [%d, %d]" % (istart_d, iend_d,
                                                    istart_s, iend_s))
    return istart_d, iend_d, istart_s, iend_s
