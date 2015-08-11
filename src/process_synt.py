import time
import glob
import os
from proc_util import process_synt_file
from source import CMTSource

syntdir = "./synt"
period_band = [27., 60.]
cmtfile = "./cmt"
outputdir = "./synt_proc"

# read in cmt
cmtsource = CMTSource.from_CMTSOLUTION_file(cmtfile)
event_time = cmtsource.cmt_time
print "Event cmt time:", event_time

# interpolation parameter
starttime = event_time
endtime = event_time + 3600.0  # 1 hour = 3600 sec
interp_deltat = 0.5 # deltat = 0.1s

synt_filelist = glob.glob(os.path.join(syntdir, "*.sac"))
print "Total number of synt files: %d" % len(synt_filelist)

# clock
t1 = time.time()

for synt_file in synt_filelist:
    process_synt_file(synt_file, period_band=period_band, interp_deltat=interp_deltat,
                      interp_starttime=interp_starttime, interp_endtime=interp_endtime,
                      outputdir=outputdir, output_format="sac", print_mode=True)

# clock
t2 = time.time()
print "Total time:", t2-t1
