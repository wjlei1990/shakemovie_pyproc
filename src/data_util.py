from obspy import read, read_inventory, readEvents
from obspy.core.util.geodetics import gps2DistAzimuth
import time
import glob
import os
import re
import numpy as np


def check_array_order(array, order="ascending"):
    """
    Check whether a array is in ascending order or descending order

    :param array:
    :return:
    """
    array = np.array(array)
    if order == "descending":
        array *= -1.0

    if array == array.sort():
        return True
    else:
        return False


def process_obsd_file(filename, stationxml_dir=None, period_band=None,
                      interp_sampling_rate=None, interp_starttime=None, interp_endtime=None,
                      outputdir=None, print_mode=False):
    """
    Processing observed data file

    :param filename:
    :param stationxml_dir:
    :param period_band:
    :param interp_sampling_rate:
    :param interp_starttime:
    :param interp_endtime:
    :param outputdir:
    :param print_mode:
    :return:
    """

    if os.path.exists(filename):
        print "file does not exist so skipped: %s" %filename
    st = read(filename)

    # calculate the frequency band based on period band(4 corners)
    period_band = np.array(period_band)
    if len(period_band) == 2:
        # if not specify the 4 corners, then use default schema
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

    print ">>>>>> Processing file: %s" % filename if print_mode else ""

    stname = st[0].stats.station
    nwname = st[0].stats.network
    # fetch station information and read stationxml file
    stationxml_file = os.path.join(stationxml_dir, "%s.%s.xml" %(nwname, stname))
    print "\tAttaching stationxml file:", stationxml_file
    inv = read_inventory(stationxml_file)
    st.attach_response(inv)

    for i, tr in enumerate(st):
        tr.detrend("linear")
        tr.detrend("demean")
        tr.taper(max_percentage=0.05, type="hann")

        tr.remove_response(output="DISP", pre_filt=pre_filt, zero_mean=False,
                           taper=False)

        tr.detrend("linear")
        tr.detrend("demean")
        tr.taper(max_percentage=0.05, type="hann")

        npts = int((interp_endtime - interp_starttime) / interp_sampling_rate)
        try:
            tr.interpolate(sampling_rate=interp_sampling_rate, starttime=interp_starttime,
                           npts=npts)
        except Exception as e:
            print "Can not interpolate data: %s" % filename
            print "Error message:", e
            return

    return st
