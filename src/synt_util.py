from obspy import read
import os
from obspy.signal.invsim import c_sac_taper
import numpy as np

def check_array_order(array, order="ascending"):
    """
    Check whether a array is in ascending order or descending order
    :param array:
    :return:
    """
    array = np.array(array)
    if order == "descending":
        array = array * -1.0

    if array == array.sort():
        return True
    else:
        return False

def filter_synt(tr, pre_filt):
    """
    Perform a frequency domain taper like during the response removal
    just without an actual response...

    :param tr:
    :param pre_filt:
    :return:
    """

    data = tr.data.astype(np.float64)

    # smart calculation of nfft dodging large primes
    from obspy.signal.util import _npts2nfft
    nfft = _npts2nfft(len(data))

    fy = 1.0 / (tr.stats.delta * 2.0)
    freqs = np.linspace(0, fy, nfft // 2 + 1)

    # Transform data to Frequency domain
    data = np.fft.rfft(data, n=nfft)
    data *= c_sac_taper(freqs, flimit=pre_filt)
    data[-1] = abs(data[-1]) + 0.0j
    # transform data back into the time domain
    data = np.fft.irfft(data)[0:len(data)]
    # assign processed data and store processing information
    tr.data = data


def process_synt_file(synt_file, period_band=None, interp_sampling_rate=None,
                      interp_starttime=None, interp_endtime=None, outputdir=None):
    """
    Processing synthetic file(without rotation)
    Operations include: 1) demean, detrend, taper; 2) filtering; 3) interpolation

    :param synt_file: synt seismogram file path
    :param period_band: period band used for filtering. Dimension could be either 2 or 4. If 2, then 2 inner corners
        would be specified; if 4, specify all 4 corners. It should be in ascending order.
    :param interp_sampling_rate:
    :param interp_starttime:
    :param interp_endtime:
    :param outputdir:
    :return:
    """

    # read in synt file
    if os.path.exists(synt_file):
        print "Synt file not exists: %s" % synt_file
    st = read(synt_file)

    # calculate the frequency band based on period band(4 corners)
    period_band = np.array(period_band)
    if len(period_band) == 2:
        freq_c2 = 1.0 / max(period_band)
        freq_c3 = 1.0 / min(period_band)
        freq_c1 = 0.80 * freq_c2
        freq_c4 = 1.25 * freq_c3
    elif len(period_band) == 4:
        if not check_array_order(period_band, order='ascending'):
            raise ValueError("Period band assigned not in ascending order: %s" %period_band)
        freq_c1 = 1.0 / period_band[3]
        freq_c2 = 1.0 / period_band[2]
        freq_c3 = 1.0 / period_band[1]
        freq_c4 = 1.0 / period_band[0]
    else:
        raise ValueError("Length of period_band should be 2 or 4")
    pre_filt = np.array([freq_c1, freq_c2, freq_c3, freq_c4])

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    nrecords = len(st)
    # processing
    for i, tr in enumerate(st):
        # get component name
        print "Processing %d of total %d" %(i+1, nrecords)
        tr.detrend("linear")
        tr.detrend("demean")
        tr.taper(max_percentage=0.05, type="hann")

        # Perform a frequency domain taper like during the response removal
        # just without an actual response...
        filter_synt(tr, pre_filt)

        tr.detrend("linear")
        tr.detrend("demean")
        tr.taper(max_percentage=0.05, type="hann")

        # interpolation
        npts = (interp_endtime - interp_starttime) / interp_sampling_rate
        try:
            tr.interpolate(sampling_rate=interp_sampling_rate, starttime=interp_starttime,
                           npts=npts)
        except Exception as e:
            print "Error: Can not interp synt"
            print "Error message: %s" % e
            return
