#name:          SonarBeamCorrection.py
#created:       May 2016
#by:            p.kennedy@fugro.com
#description:   This script will compute a waterfall image of the raw sonar data into a greyscale png file.  Very nice for a quick QC of data quality
#notes:         See main at end of script for example how to use this
#based on XTF version 26 18/12/2008
#version 1.00

#2DO
# ensure the waterfall works when files without navigation are present in the file
# add support for altitude based BC tables.  This means compute the altitude range, then segment on a regular intervel (default 1m)
# add support for user selection of palette.  Use the edgetech discover palettes files.  default to be greyscale
# add support for logarithmic color ramp
# add support for color ramp scaling to the command line (min/max value)
# add support for channel selection, 0,1 or 2,3  default is 0,1
# what to do if altitude in XTF == 0.  needs testing for div zero.
# add support for slant range correction, ie remove the water column from the resulting waterfall images

#DONE
# revert to local copy of pyxtf.py so tool is more portable
# designed to handle both geographic and grid positions seamlessly
# creates a nadir gap polygon representing the null area under the sidescan which is poor quality data.
# uses the sensor position, altitude and rather than computing CMG heading which is a bit wobbly when navigation data has duplicate points
# add support for stretch and decimation instead of auto computing based on navigation.  We may well not have navigation!
# -w option now produces a nice png file of the sonar
# now uses the pyxtf project rather than a local copy
# moved over to revised XTFreader which reads packets rather than records and handles padbytes
# improved user feedback
# performance improvements
# current performance is 27 seconds to process 8000 pings (half an hour of AUV collection.  This is ok
# added a simple waterfall image creation funciton.  We can use this to confirm the beam correction is valid
# auto compute the approximate Y scale form the sample rate so the waterfall image is approximately isometric
# do we need to take into account changes in sample frequency.  This will screw up the BC table 
# based on XTF version 26 18/12/2008
# based on existing SonarCoverage.py script
# initial implementation
=======

