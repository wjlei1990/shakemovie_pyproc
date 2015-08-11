from obspy import readEvents
import time
import glob
import os
import sys
from source import CMTSource
from proc_util import process_obsd_file

datadir = "../testdata/obsd"
period_band = [27., 60.]
cmtfile = "../testdata/cmt/C201311120703A"
outputdir = "../testdata/obsd_proc"
stationxmldir = "../testdata/stationxml"

# read in cmt
cmtsource = CMTSource.from_CMTSOLUTION_file(cmtfile)
event_time = cmtsource.cmt_time
print "Event cmt time:", event_time

# interpolation parameter
starttime = event_time
endtime = event_time + 3600.0  # 1 hour = 3600 sec
interp_deltat = 0.5  # deltat = 0.1s

obsd_filelist = glob.glob(os.path.join(datadir, "*.mseed"))
print "Total number of observed seismograms: %d" % len(obsd_filelist)

# clock
t1 = time.time()

for obsd_file in obsd_filelist:
    process_obsd_file(obsd_file, stationxmldir=stationxmldir, period_band=period_band,
                      interp_deltat=interp_deltat, interp_starttime=starttime, interp_endtime=endtime,
                      outputdir=outputdir, output_format="SAC", print_mode=True)

# clock
t2 = time.time()
print "Total time:", t2-t1
