import obspy
from source import CMTSource
from proc_util import rotate_trace

inputdir = "../testdata/synt_proc"
channel_code="MX"
outputdir = "../testdata/rotate"
cmtfile = "../testdata/cmt/C201311120703A"
staxmldir = "../testdata/stationxml"

cmtsource = CMTSource.from_CMTSOLUTION_file(cmtfile)

rotate_trace(inputdir, data_format="sac", channel_code=channel_code, outputdir=outputdir, event=cmtsource, stationxmldir=staxmldir)

#def rotate_seismogram(inputdir, data_format="sac", channel_code="MX", outputdir=None, outputformat="sac",
#                      event=None, stationxmldir=None):
