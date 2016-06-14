#name:          SonarBeamCorrection.py
#created:       May 2016
#by:            p.kennedy@fugro.com
#description:   This script will compute a beam correction file from a series of XTF files and save the results to a CSV file
#notes:         See main at end of script for example how to use this
#version 1.00

#2DO
# auto compute the approximate Y scale form the sample rate so the waterfall image is approximately isometric
# add support for altitude based BC tables.  This means compute the altitude range, then segment on a regular intervel (default 1m)

#DONE
# do we need to take into account changes in sample frequency.  THis will screw up the BC table 
# based on XTF version 26 18/12/2008
# based on existing SonarCoverage.py script
# initial implementation

import math
import argparse
import sys
import csv
import pyXTF
import geodetic
import os
import datetime
from glob import glob
import time
import numpy as np
from PIL import Image

# we calculate teh angle in degrees from the altitude and sample number
def calcAngleFromSample(altitude):
    return (altitude * 0.70) / 2.0

def calculateRangeBearingFromPosition(easting1, northing1, easting2, northing2):
    """given 2 east, north, pairs, compute the range and bearing"""

    dx = easting2-easting1
    dy = northing2-northing1

    bearing = 90 - (180/math.pi)*math.atan2(northing2-northing1, easting2-easting1)
    return (math.sqrt((dx*dx)+(dy*dy)), bearing)

def main():

    start_time = time.time() # time the process
    parser = argparse.ArgumentParser(description='Read XTF file and create either a coverage or Nadir gap polygon.')
    parser.add_argument('-c', action='store_true', default=False, dest='createBC', help='-c compute a new Beam Correction file based on contents of XTF file[s]')
    parser.add_argument('-i', dest='inputFile', action='store', help='-i <XTFfilename> input XTF filename to analyse')
    parser.add_argument('-w', dest='createWaterfall', default = False, action='store_true', help='-w create a waterfall image from the XTF file')
    parser.add_argument('-o', dest='outputFile', action='store', help='-o <filename> output Beam Correction filename.  [default = BC.csv]')
    parser.add_argument('-odix', dest='outputFolder', action='store', help='-odix <folder> output folder to store Beam Correction file.  If not specified, the files will be alongside the input XTF file')
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
   
    if args.outputFolder == None:
        firstFile = glob(args.inputFile)[0]
        args.outputFolder = os.path.abspath(os.path.join(firstFile, os.pardir))


    leftSide = [] #storage for the left transducer beam correction
    rightSide = [] #storage for the right transducer beam correction
    maxSamplesPort = 0
    maxSamplesStbd = 0
    minAltitude = 99999
    maxAltitude = 0
    segmentInterval = 5 #the altitude segmentaiton interval
    
    print ("files to Process:", glob(args.inputFile))
    print ("iterating through all input files to compute the maximm size of the beam correction table.  This takes into account changes in range")
    if args.createBC:       
        for filename in glob(args.inputFile):
            samplesPort, samplesStbd, minAltitude, maxAltitude, pingCount = getSampleRange(filename)
            maxSamplesPort = max(maxSamplesPort, samplesPort)
            maxSamplesStbd = max(maxSamplesStbd, samplesStbd)
        print ("maxSamplesPort %s maxSamplesStbd %s minAltitude %s maxAltitude %s" % (maxSamplesPort, maxSamplesStbd, minAltitude, maxAltitude))
        
        numSegments = int(maxAltitude / segmentInterval) + 1 # need to add extra for zero index
        samplesPortSum = np.zeros((numSegments, maxSamplesPort), dtype=np.int32  )
        samplesPortCount = np.zeros((numSegments,maxSamplesPort), dtype=np.int  )
        samplesStbdSum = np.zeros((numSegments, maxSamplesStbd), dtype=np.int32  )
        samplesStbdCount = np.zeros((numSegments, maxSamplesStbd), dtype=np.int  )

    for filename in glob(args.inputFile):
        if args.createBC:       
            samplesPortSum, samplesPortCount, samplesStbdSum, samplesStbdCount = computeBC(filename, samplesPortSum, samplesPortCount, samplesStbdSum, samplesStbdCount, segmentInterval)
        if args.createWaterfall:       
            createWaterfall(filename)
        else:
            print ("option not yet implemented!.  Try '-n' to compute nadir gaps")
            exit (0)

    if args.createBC:       
        # now save the results to a csv file
        if args.outputFile is None:
            baseName = os.path.basename(os.path.splitext(glob(args.inputFile)[0])[0])
            if args.outputFolder is None:
                args.outputFolder = os.path.curdir(glob(args.inputFile)[0])
        else:
            baseName = os.path.basename(os.path.splitext(glob(args.outputFile)[0]))
            if args.outputFolder is None:
                args.outputFolder = os.path.curdir(glob(args.outputFile)[0])


        BCPortFile = os.path.join(args.outputFolder, baseName + "_port.csv")
        BCStbdFile = os.path.join(args.outputFolder, baseName + "_stbd.csv")

        print("saving Port CSV file:", BCPortFile)
        samplesPortAVG = np.divide(samplesPortSum, samplesPortCount)
        np.savetxt(BCPortFile, samplesPortAVG, delimiter=',')

        print("saving Stbd CSV file:", BCStbdFile)
        samplesStbdAVG = np.divide(samplesStbdSum, samplesStbdCount)
        np.savetxt(BCStbdFile, samplesStbdAVG, delimiter=',')

    print("--- %s seconds ---" % (time.time() - start_time)) # print the processing time.

    return (0)

def mergeImages(image1, image2):
    """Merge two images into one, displayed side by side
    :param file1: path to first image file
    :param file2: path to second image file
    :return: the merged Image object
    """

    (width1, height1) = image1.size
    (width2, height2) = image2.size

    result_width = width1 + width2
    result_height = max(height1, height2)

    result = Image.new('RGB', (result_width, result_height))
    result.paste(im=image1, box=(0, 0))
    result.paste(im=image2, box=(width1, 0))
    return result

def getSampleRange(filename):
    maxSamplesPort = 0
    maxSamplesStbd = 0
    minAltitude = 99999
    maxRange = 0
    maxAltitude = 0
    pingCount = 0
    #   open the XTF file for reading 
    r = pyXTF.XTFReader(filename)
    print ("computing range for file:", filename)
    while r.moreData():
        pingHdr = r.readPing()
        maxSamplesPort = max(pingHdr.pingChannel[0].NumSamples, maxSamplesPort)
        maxSamplesStbd = max(pingHdr.pingChannel[1].NumSamples, maxSamplesStbd)
        minAltitude = min(minAltitude, pingHdr.SensorPrimaryAltitude)
        maxAltitude = max(maxAltitude, pingHdr.SensorPrimaryAltitude)
        maxRange = max(maxRange, pingHdr.pingChannel[0].SlantRange)
        
        d = datetime.date(pingHdr.Year, pingHdr.Month, pingHdr.day)
        t = datetime.time(pingHdr.Hour, pingHdr.Minute, pingHdr.Second)
        dt = datetime.datetime.combine(d, t)
        print 'dt:', dt

        pingCount = pingCount + 1
    return maxSamplesPort, maxSamplesStbd, minAltitude, maxAltitude, maxRange, pingCount

# compute the segment from the altitude in a standard manner.  
def altitudeToSegment(altitude, segmentInterval):    
    segmentnumber = int( altitude / segmentInterval)
    # segmentnumber = round(altitude - segmentInterval/2,-1) / segmentInterval
    return segmentnumber

# create a simple waterfall image from the sonar data using standard nunmpy and PIL python pachakages
def createWaterfall(filename):
    yStretch = 2
    decimation = 4
    #compute the size of the array we will need to create
    maxSamplesPort, maxSamplesStbd, minAltitude, maxAltitude, maxSlantRange, pingCount= getSampleRange(filename)
    
    acrossTrackSampleInterval = (maxSlantRange / maxSamplesPort) * decimation # sample interval in metres
    
    alongTrackSampleInterval = 2 # assume 4 knots
    
    yStretch = math.ceil(alongTrackSampleInterval / acrossTrackSampleInterval)
    
    # portChannel = np.zeros((pingCount, maxSamplesPort), dtype=np.int16  )
    # stbdChannel = np.zeros((pingCount, maxSamplesStbd), dtype=np.int16  )

    #   open the XTF file for reading 
    print ("Opening file:", filename)
    r = pyXTF.XTFReader(filename)
    # ping = 0
    pc = []
    sc = []
    
    print ("Loading Data...")
    while r.moreData():
        pingHdr = r.readPing()
        for i in range (yStretch):
            pc.append(pingHdr.pingChannel[0].data[::decimation])            
            sc.append(pingHdr.pingChannel[1].data[::decimation])            
        # ping = ping + 1

    print ("Converting to Image...")
    portImage = samplesToImage(pc)
    stbdImage = samplesToImage(sc)

    #merge the images into a single image and save to disc
    print ("Merging Channels to single image...")
    mergedImage = mergeImages(portImage, stbdImage).convert('L')
    mergedImage.show()
    print ("Saving Image...")    
    mergedImage.save(os.path.splitext(filename)[0]+'.png')
        
def samplesToImage(samples):
    zg_LL = 5 # min and max grey scales
    zg_UL = 250

    #create numpy arrays so we can compute stats
    channel = np.array(samples)
    
    # now filter for outliers
    # import scipy
    
    channelMin = channel.min()
    channelMax = channel.max()
    
    if channelMin > 0:
        zs_LL = math.log(channelMin)
    else:
        zs_LL = 0
    if channelMax > 0:
        zs_UL = math.log(channelMax)
    else:
        zs_UL = 0
    
    if (zs_UL - zs_LL) is not 0:
        conv_01_99 = ( zg_UL - zg_LL ) / ( zs_UL - zs_LL )
   
    #we can expect some divide by zero errors, so suppress 
    np.seterr(divide='ignore')
    channel = np.log(samples)
    channel = np.subtract(channel, zs_LL)
    channel = np.multiply(channel, conv_01_99)
    channel = np.subtract(zg_UL, channel)
    ch = channel.astype('uint8')
    image = Image.fromarray(ch).convert('L')
    
    from PIL import ImageFilter
    im2 = image.filter(ImageFilter.MedianFilter(3))

    return im2
    
    
    
# maxSamples lets us know how wide to make the beam correciton table.  This is computed by a pre-iteration stage
def computeBC(filename, samplesPortSum, samplesPortCount, samplesStbdSum, samplesStbdCount, segmentInterval):

    #   open the XTF file for reading 
    print ("Opening file:", filename)
    r = pyXTF.XTFReader(filename)
    while r.moreData():
        pingHdr = r.readPing()
        segment = altitudeToSegment(pingHdr.SensorPrimaryAltitude, segmentInterval)
        
        samplesPortSum[segment] = np.add(samplesPortSum[segment], pingHdr.pingChannel[0].data)            
        samplesPortCount[segment] = np.add(samplesPortCount[segment],1)            
        
        samplesStbdSum[segment] = np.add(samplesStbdSum[segment], pingHdr.pingChannel[1].data)            
        samplesStbdCount[segment] = np.add(samplesStbdCount[segment],1)            

        if pingHdr.PingNumber % 1000 == 0:
            print ("Ping: %f" % (pingHdr.PingNumber))                  
    
    print("Complete reading XTF file :-)")
    return (samplesPortSum, samplesPortCount, samplesStbdSum, samplesStbdCount)
    
if __name__ == "__main__":
    main()
