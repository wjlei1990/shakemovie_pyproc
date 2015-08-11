from obspy import read, read_inventory
import os
from obspy.signal.invsim import c_sac_taper
import numpy as np
from obspy.signal.util import _npts2nfft


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


def write_out_seismogram(st, outputdir, output_format, originfn=None):
    output_format = output_format.lower()
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if output_format == "mseed":
        outputfn = os.path.basename(originfn)
        outputpath = os.path.join(outputdir, outputfn+".proc")
        st.write(outputpath, format="mseed")
    if output_format == "sac":
        for tr in st:
            station = tr.stats.station
            network = tr.stats.network
            channel = tr.stats.channel
            location = tr.stats.location
            outputfn = "%s.%s.%s.%s.sac" % (network, station, location, channel)
            outputpath = os.path.join(outputdir, outputfn)
            tr.write(outputpath, format="sac")


def process_synt_file(filename, period_band=None, interp_deltat=None,
                      interp_starttime=None, interp_endtime=None, outputdir=None,
                      output_format="SAC", print_mode=False):
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
    if not os.path.exists(filename):
        print "file does not exist so skipped: %s" % filename
        return
    st = read(filename)

    # calculate the frequency band based on period band(4 corners)
    period_band = np.array(period_band)
    if len(period_band) == 2:
        # if not specify the 4 corners, then use default schema
        freq_c2 = 1.0 / max(period_band)
        freq_c3 = 1.0 / min(period_band)
        freq_c1 = 0.8 * freq_c2
        freq_c4 = 1.2 * freq_c3
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

    nrecords = len(st)
    # processing
    print ">>>>>> Processing file: %s" % filename if print_mode else ""
    for i, tr in enumerate(st):
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
        npts = int((interp_endtime - interp_starttime) / interp_deltat)
        try:
            tr.interpolate(sampling_rate=1.0/interp_deltat, starttime=interp_starttime,
                           npts=npts)
        except Exception as e:
            print "Error: Can not interp synt"
            print "Error message: %s" % e
            return

        try:
            tr.data /= tr.stats.calib
        except:
            pass

    if outputdir is not None:
        write_out_seismogram(st, outputdir, output_format)

    return st


def process_obsd_file(filename, stationxmldir=None, period_band=None,
                      interp_deltat=None, interp_starttime=None, interp_endtime=None,
                      outputdir=None, output_format="SAC", print_mode=False):
    """
    Processing observed data file

    :param filename:
    :param stationxmldir:
    :param period_band:
    :param interp_sampling_rate:
    :param interp_starttime:
    :param interp_endtime:
    :param outputdir:
    :param print_mode:
    :return:
    """

    if not os.path.exists(filename):
        print "file does not exist so skipped: %s" % filename
        return
    st = read(filename)

    # calculate the frequency band based on period band(4 corners)
    period_band = np.array(period_band)
    if len(period_band) == 2:
        # if not specify the 4 corners, then use default schema
        freq_c2 = 1.0 / max(period_band)
        freq_c3 = 1.0 / min(period_band)
        freq_c1 = 0.8 * freq_c2
        freq_c4 = 1.2 * freq_c3
    elif len(period_band) == 4:
        if not check_array_order(period_band, order='ascending'):
            raise ValueError("Period band assigned not in ascending order: %s" % period_band)
        freq_c1 = 1.0 / period_band[3]
        freq_c2 = 1.0 / period_band[2]
        freq_c3 = 1.0 / period_band[1]
        freq_c4 = 1.0 / period_band[0]
    else:
        raise ValueError("Length of period_band should be 2 or 4")
    #pre_filt = np.array([freq_c1, freq_c2, freq_c3, freq_c4])
    pre_filt = (freq_c1, freq_c2, freq_c3, freq_c4)

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    print ">>>>>> Processing file: %s" % filename if print_mode else ""

    stname = st[0].stats.station
    nwname = st[0].stats.network
    # fetch station information and read stationxml file
    stationxml_file = os.path.join(stationxmldir, "%s.%s.xml" %(nwname, stname))
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

        npts = int((interp_endtime - interp_starttime) / interp_deltat)
        try:
            tr.interpolate(sampling_rate=1.0/interp_deltat, starttime=interp_starttime,
                           npts=npts)
        except Exception as e:
            print "Can not interpolate data: %s" % filename
            print "Error message:", e
            return

    if outputdir is not None:
        write_out_seismogram(st, outputdir, output_format, originfn=filename)

    return st
