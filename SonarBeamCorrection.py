#name:          SonarBeamCorrection.py
#created:       May 2016
#by:            p.kennedy@fugro.com
#description:   This script will compute a beam correction file from a series of XTF files and save the results to a CSV file
#notes:         See main at end of script for example how to use this
#version 1.00

#2DO
# add support for altitude based BC tables.  This means compute the altitude range, then segment on a regular intervel (default 1m)
# add support for user selection of palette.  Use the edgetech discover palettes files.  default to be greyscale
# add support for logarithmic color ramp
# add support for color ramp scaling to the command line (min/max value)
# add support for channel selection, 0,1 or 2,3  default is 0,1

#DONE
# add support for stretch and decimation instead of auto computing based on navigation.  We may well not have navigation!
# -w option now produces a nice png file of the sonar
# now uses the pyxtf project rather than a local copy
# moved over to revised XTFreader which reads packets rather than records and handles padbytes
# improved user feedbakc
# performance improvements
# current performane is 27 econds to process 8000 pings (half an hour of AUV collection.  This is ok
# added a simple waterfall image creation funciton.  We can use this to confirm the beam correction is valid
# auto compute the approximate Y scale form the sample rate so the waterfall image is approximately isometric
# do we need to take into account changes in sample frequency.  THis will screw up the BC table 
# based on XTF version 26 18/12/2008
# based on existing SonarCoverage.py script
# initial implementation


import bisect
import math
import argparse
import sys
# good for debugging as we only have the 1 file, but for production we should move these into the same folder as the script so it is all self contained.
sys.path.append('../pyxtf')
import csv
import pyXTF
import geodetic
import os
from datetime import datetime
from glob import glob
import time
import numpy as np
from PIL import Image
from scipy import interpolate

# from pylab import imshow, show, get_cmap
# import matplotlib
# import pandas as pd
import matplotlib.pyplot as plt

def main():
    start_time = time.time() # time the process
    parser = argparse.ArgumentParser(description='Read XTF file and create either a coverage or Nadir gap polygon.')
    parser.add_argument('-a', action='store_true', default=False, dest='applyBC', help='-a : The Beam Correction file to apply when processing the XTF file[s]. Use the -c option to create the BCFile [default = False]')
    parser.add_argument('-c', action='store_true', default=False, dest='createBC', help='-c : Compute a new Beam Correction file based on contents of XTF file[s]. File create will be called beamcorreciton_port and beamcorrection_stbd.csv')
    parser.add_argument('-i', dest='inputFile', action='store', help='-i <XTFfilename> : input XTF filename to analyse')
    parser.add_argument('-color', dest='color', default = 'graylog', action='store', help='-color <paletteName> : Specify the color palette.  Options are : -color yellow_brown_log, -color gray, -color yellow_brown or any of the palette filenames in the script folder. [Default = graylog for a grayscale logarithmic palette.]' )
    parser.add_argument('-clip', dest='clip', default = 0, action='store', help='-clip <value> : Clip the minimum and maximum edges of the data by this percentage so the color stretch better represents the data.  [Default - 0.  A good value is -clip 5.]')
    parser.add_argument('-invert', dest='invert', default = False, action='store_true', help='-invert : Inverts the color palette')
    parser.add_argument('-o', dest='outputFile', action='store', help='-o <filename> : Output Beam Correction filename.  [default = BC.csv]')
    parser.add_argument('-odix', dest='outputFolder', action='store', help='-odix <folder> : Output folder to store Beam Correction file.  If not specified, the files will be alongside the input XTF file')
    parser.add_argument('-w', dest='createWaterfall', default = False, action='store_true', help='-w : Create a waterfall image from the XTF file')
    parser.add_argument('-x', action='store', default=4, dest='decimation', help='-x <value> : Decimate the samples in the across track direction.  For full resolution, use 1.  For a faster process, use 10. The process will look at the data and inform you what the X axis pixel size represents with this decimation factor. [Default = 4.]')
    parser.add_argument('-y', action='store', default=-1, dest='stretch', help='-y <value> : Set the alongtrack scale factor. The process will look at the data and inform you what the Y axis pixel size repesents with this y scale. [Default = -1 for auto stretch.]')
    parser.add_argument('-lf', action='store_true', default=True, dest='processLFData', help='-lf : process the LOW frequency data in the file (channels 0 & 1. [Default = True.]')
    parser.add_argument('-hf', action='store_true', default=False, dest='processHFData', help='-hf : process the HIGH frequency data in the file (channels 2 & 3. [Default = False.]')
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()

    print ("processing with settings: ", args)

    # default to process the low frequency channels which are 0 an 1
    channelA = 0
    channelB = 1        
    if args.processHFData:
        channelA = 2
        channelB = 3
                
    if args.outputFolder == None:
        firstFile = glob(args.inputFile)[0]
        args.outputFolder = os.path.abspath(os.path.join(firstFile, os.pardir))

    leftSide = [] #storage for the left transducer beam correction
    rightSide = [] #storage for the right transducer beam correction
    maxSamplesPort = 0
    maxSamplesStbd = 0
    minAltitude = 99999
    maxAltitude = 0
    segmentInterval = 100 #the altitude segmentaiton interval
    
    print ("Files to Process:", glob(args.inputFile))
    print ("Iterating through all input files to compute the maximum size of the beam correction table.  This takes into account changes in range")
    if args.createBC:       
        for filename in glob(args.inputFile):
            samplesPort, samplesStbd, minAltitude, maxAltitude, maxRange, pingCount, meanSpeed, navigation = getSampleRange(filename, channelA, channelB, False)
            maxSamplesPort = max(maxSamplesPort, samplesPort)
            maxSamplesStbd = max(maxSamplesStbd, samplesStbd)
        print ("maxSamplesPort %s maxSamplesStbd %s minAltitude %s maxAltitude %s" % (maxSamplesPort, maxSamplesStbd, minAltitude, maxAltitude))
        
        numSegments = int(maxAltitude / segmentInterval) + 1 # need to add extra for zero index
        samplesPortSum = np.zeros((numSegments, 90), dtype=np.int32  )
        samplesPortCount = np.zeros((numSegments,90), dtype=np.int  )
        samplesStbdSum = np.zeros((numSegments, 90), dtype=np.int32  )
        samplesStbdCount = np.zeros((numSegments, 90), dtype=np.int  )

        # samplesPortSum = np.zeros((numSegments, maxSamplesPort), dtype=np.int32  )
        # samplesPortCount = np.zeros((numSegments,maxSamplesPort), dtype=np.int  )
        # samplesStbdSum = np.zeros((numSegments, maxSamplesStbd), dtype=np.int32  )
        # samplesStbdCount = np.zeros((numSegments, maxSamplesStbd), dtype=np.int  )

    BCPortFile = ""
    BCStbdFile = ""
    if args.applyBC:
        folder = os.path.dirname(glob(args.inputFile)[0])
        baseName = "beamcorrection"
        BCPortFile = os.path.join(folder, baseName + "_port.csv")
        BCStbdFile = os.path.join(folder, baseName + "_stbd.csv")

    for filename in glob(args.inputFile):
        if args.createBC:
            samplesPortSum, samplesPortCount, samplesStbdSum, samplesStbdCount = computeBC(filename, channelA, channelB, samplesPortSum, samplesPortCount, samplesStbdSum, samplesStbdCount, segmentInterval)
        if args.createWaterfall:       
            createWaterfall(filename, channelA, channelB, args.invert, args.color, float(args.clip), int(args.decimation), int(args.stretch), int(args.applyBC), BCPortFile, BCStbdFile, segmentInterval)

    if args.createBC:       
        # now save the results to a csv file
        if args.outputFile is None:
            baseName = os.path.basename(os.path.splitext(glob(args.inputFile)[0])[0])
            if args.outputFolder is None:
                args.outputFolder = os.path.dirname(glob(args.inputFile)[0])
        else:
            baseName = os.path.basename(os.path.splitext(glob(args.outputFile)[0]))
            if args.outputFolder is None:
                args.outputFolder = os.path.dirname(glob(args.outputFile)[0])

        # BCPortFile = os.path.join(args.outputFolder, baseName + "_port.csv")
        # BCStbdFile = os.path.join(args.outputFolder, baseName + "_stbd.csv")
        baseName = "beamcorrection"
        BCPortFile = os.path.join(args.outputFolder, baseName + "_port.csv")
        BCStbdFile = os.path.join(args.outputFolder, baseName + "_stbd.csv")

        print("Saving Port CSV file:", BCPortFile)
        samplesPortAVG = np.divide(samplesPortSum, samplesPortCount)
        # apply a smoothing fileter to the data.  There is no value having spikes in the beam correction file!
        # samplesPortAVG = np.apply_along_axis(geodetic.medfilt, 1, samplesPortAVG, 9)
        np.savetxt(BCPortFile, samplesPortAVG, fmt='%.3f', delimiter=',')

        print("Saving Stbd CSV file:", BCStbdFile)
        samplesStbdAVG = np.divide(samplesStbdSum, samplesStbdCount)
        # samplesStbdAVG = np.apply_along_axis(geodetic.medfilt, 1, samplesStbdAVG, 9)
        np.savetxt(BCStbdFile, samplesStbdAVG, fmt='%.3f', delimiter=',')

    print("Process complete. Duration: %.2f s" % (time.time() - start_time)) # print the processing time.

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

    result = Image.new('L', (result_width, result_height))
    result.paste(im=image1, box=(0, 0))
    result.paste(im=image2, box=(width1, 0))
    return result

def getSampleRange(filename, channelA, channelB, loadNavigation):
    """iterate through the file to find the extents for range, time and samples.  These are all needed in subsequent processing """
    maxSamplesPort = 0
    maxSamplesStbd = 0
    minAltitude = 99999
    maxRange = 0
    maxAltitude = 0
    pingCount = 0
    pingTimes = []
    navigation = 0
    
    print("Gathering data limits...")
    #   open the XTF file for reading 
    r = pyXTF.XTFReader(filename)
    if loadNavigation:
        navigation = r.loadNavigation()
    # meanSpeed, navigation = r.computeSpeedFromPositions(navigation)
    meanSpeed = 1
    start_time = time.time() # time the process
    
    while r.moreData():
        ping = r.readPacket()
        maxSamplesPort = max(ping.pingChannel[channelA].NumSamples, maxSamplesPort)
        maxSamplesStbd = max(ping.pingChannel[channelB].NumSamples, maxSamplesStbd)
        minAltitude = min(minAltitude, ping.SensorPrimaryAltitude)
        maxAltitude = max(maxAltitude, ping.SensorPrimaryAltitude)
        maxRange = max(maxRange, ping.pingChannel[channelA].SlantRange)
        pingCount = pingCount + 1

    print("Get Sample Range Duration %.3fs" % (time.time() - start_time)) # print the processing time.
    return maxSamplesPort, maxSamplesStbd, minAltitude, maxAltitude, maxRange, pingCount, meanSpeed, navigation

def altitudeToSegment(altitude, segmentInterval):    
    """Compute the segment from the altitude in a standard manner"""  
    segmentnumber = int( altitude / segmentInterval)
    # segmentnumber = round(altitude - segmentInterval/2,-1) / segmentInterval
    return segmentnumber
    
##############################################################################################################################
def createWaterfall(filename, channelA, channelB, invert, colorScale, clip, decimation, stretch, applyBC, BCPortFile, BCStbdFile, segmentInterval):
    """create a simple waterfall image from the sonar data using standard nunmpy and PIL python pachakages"""
    
    if applyBC:
        if os.path.isfile(BCPortFile):
            PortBC = np.loadtxt(BCPortFile, delimiter=',', dtype=float)
            # take every nth sample using the specified decimation factor
            # PortBC = PortBC[::1,::decimation]
            # PortBC = PortBC[::decimation]
        else:
            applyBC = False
        if os.path.isfile(BCStbdFile):
            StbdBC = np.loadtxt(BCStbdFile, delimiter=',')
            # StbdBC = StbdBC[::1,::decimation]
            # StbdBC = StbdBC[::decimation]
        else:
            applyBC = False
            
    #compute the size of the array we will need to create
    maxSamplesPort, maxSamplesStbd, minAltitude, maxAltitude, maxSlantRange, pingCount, meanSpeed, navigation = getSampleRange(filename, channelA, channelB, True)
    
    acrossTrackSampleInterval = (maxSlantRange / maxSamplesPort) * decimation # sample interval in metres
    
    # to make the image somewhat isometric, we need to compute the alongtrack sample interval.  this is based on the ping times, number of pings and mean speed  where distance = speed * duration
    distance = meanSpeed * (navigation[-1].dateTime.timestamp() - navigation[0].dateTime.timestamp())
    alongTrackSampleInterval = (distance / pingCount) 
    
    stretch = math.ceil(alongTrackSampleInterval / acrossTrackSampleInterval)
    
    print ("Creating waterfall with pixel resolution acrosstrack %.3f alongtrack %.3f" % (acrossTrackSampleInterval, alongTrackSampleInterval/stretch))
    
    #   open the XTF file for reading 
    print ("Opening file:", filename)
    r = pyXTF.XTFReader(filename)
    pc = []
    sc = []
    
    print ("Loading Data...")
    
    previousRange = 0
    while r.moreData():
        ping = r.readPacket()
        # this is not a ping so skip it
        if ping == -999:
            continue
        
        # Using the angular correction file, calculate a per pixel beam correction for this range / numsamples.  We can reuse this until the range or numsamples changes. This is far more efficient than computing a take off angle for every sample in every ping
        if applyBC:
            if not (previousRange == int (ping.pingChannel[0].SlantRange)):
                previousRange = int (ping.pingChannel[0].SlantRange) 
                PortARC = []
                StbdARC = []
                print ("Computing angular response from the Beam Correction tables for slantrange:", ping.pingChannel[0].SlantRange)
                for ARC in PortBC:
                    PortARC.append(angularResponseToPing(ARC, ping.pingChannel[0].SlantRange, int (ping.pingChannel[0].NumSamples / decimation), ping.SensorPrimaryAltitude))
                portLift = np.average(PortARC)
        
        # now do the port channel
        channel = np.array(ping.pingChannel[0].data[::decimation])
        channel = np.multiply(channel, math.pow(2, -ping.pingChannel[0].Weight))
        rawPortData = channel.tolist()
        channelCorrected = channel
        channelCorrected = geodetic.medfilt(channelCorrected, 5)
        # apply the beam correction here
        if applyBC:
            segment = altitudeToSegment(ping.SensorPrimaryAltitude, segmentInterval)
            # correction = PortBC
            correction = PortARC[segment]
            channelCorrected = np.subtract(channel, correction)
            # channelCorrected = np.divide(channel, correction)
            # channelCorrected = np.add(channel, portLift)            
            filteredPortData = channelCorrected.tolist()

            plt.plot(channel, label="Raw Channel")
            plt.plot(channelCorrected, label="corrected")
            plt.plot(correction, label="correction", color='red')
            # # plt.plot(PortBC[segment-1], label="Beam correction -1")
            plt.xlabel('sample')
            plt.ylabel('intensity')
            plt.grid(True)
            # plt.subplots.ax.set_color_cycle(['red', 'black', 'yellow'])
            plt.legend()
            plt.show()            

        # now do the stbd channel
        channel = np.array(ping.pingChannel[1].data[::decimation])
        channel = np.multiply(channel, math.pow(2, -ping.pingChannel[1].Weight))
        rawStbdData = channel.tolist()
        channelCorrected = channel
        channelCorrected = geodetic.medfilt(channelCorrected, 5)
        # if applyBC:
        #     segment = altitudeToSegment(ping.SensorPrimaryAltitude, segmentInterval)
        #     channelCorrected = np.subtract(channel, StbdBC[segment])            
        # filteredStbdData = channelCorrected.tolist()

        # Add the data to the waterfall, stretching along the Y axis as we go so the image is close to isometric 
        for i in range (stretch):
            # pc.insert(0, filteredPortData)
            pc.insert(0, rawPortData)
            sc.insert(0, rawStbdData)

    portImage = []
    stbdImage = []
    
    if colorScale.lower() == "graylog": 
        print ("Converting to Image with graylog scale...")
        portImage = samplesToGrayImageLogarithmic(pc, invert, clip)
        stbdImage = samplesToGrayImageLogarithmic(sc, invert, clip)
    elif colorScale.lower() == "gray":
        print ("Converting to Image with gray scale...")
        portImage = samplesToGrayImage(pc, invert, clip)
        stbdImage = samplesToGrayImage(sc, invert, clip)
    else:
        print ("Converting to Image with default gray log scale...")
        portImage = samplesToColorImage(pc, invert, clip, colorScale)
        stbdImage = samplesToColorImage(sc, invert, clip, colorScale)
        
    # merge the images into a single image and save to disc
    print ("Merging Channels to single image...")
    mergedImage = mergeImages(portImage, stbdImage)
    print ("Saving Image...")    
    mergedImage.save(os.path.splitext(filename)[0]+'.png')
    print ("Image saved to:", os.path.splitext(filename)[0]+'.png')    
    
##################################################################################################################
def findMinMaxClipValues(channel, clip):
    print ("Clipping data with an upper and lower percentage of:", clip)
    # compute a histogram of teh data so we can auto clip the outliers
    bins = np.arange(np.floor(channel.min()),np.ceil(channel.max()))
    values, base = np.histogram(channel, bins=bins, density=1)    

    # instead of spreading across the entire data range, we can clip the outer n percent by using the cumsum.
    # from the cumsum of histogram density, we can figure out what cut off sample amplitude removes n % of data
    cumsum = np.cumsum(values)   
    
    minimumBinIndex = bisect.bisect(cumsum,clip/100)
    maximumBinIndex = bisect.bisect(cumsum,(1-clip/100))

    return minimumBinIndex, maximumBinIndex

###################################
# zg_LL = lower limit of grey scale
# zg_UL = upper limit of grey scale
# zs_LL = lower limit of samples range
# zs_UL = upper limit of sample range
def samplesToGrayImage(samples, invert, clip):
    zg_LL = 5 # min and max grey scales
    zg_UL = 250
    zs_LL = 0 
    zs_UL = 0
    conv_01_99 = 1
    
    #create numpy arrays so we can compute stats
    channel = np.array(samples)   

    # compute the clips
    if clip > 0:
        zs_LL, zs_UL = findMinMaxClipValues(channel, clip)
    else:
        zs_LL = channel.min()
        zs_UL = channel.max()
    
    # this scales from the range of image values to the range of output grey levels
    if (zs_UL - zs_LL) is not 0:
        conv_01_99 = ( zg_UL - zg_LL ) / ( zs_UL - zs_LL )
   
    #we can expect some divide by zero errors, so suppress 
    np.seterr(divide='ignore')
    # channel = np.log(samples)
    channel = np.subtract(channel, zs_LL)
    channel = np.multiply(channel, conv_01_99)
    if invert:
        channel = np.subtract(zg_UL, channel)
    else:
        channel = np.add(zg_LL, channel)
    image = Image.fromarray(channel).convert('L')
    return image

def samplesToColorImage(samples, invert, clip, palette):
    zg_LL = 5 # min and max grey scales
    zg_UL = 250
    zs_LL = 0 
    zs_UL = 0
    conv_01_99 = 1
    
    #create numpy arrays so we can compute stats
    channel = np.array(samples)   

    # compute the clips
    if clip > 0:
        zs_LL, zs_UL = findMinMaxClipValues(channel, clip)
    else:
        zs_LL = channel.min()
        zs_UL = channel.max()
    
    # this scales from the range of image values to the range of output grey levels
    if (zs_UL - zs_LL) is not 0:
        conv_01_99 = ( zg_UL - zg_LL ) / ( zs_UL - zs_LL )
   
    #we can expect some divide by zero errors, so suppress 
    np.seterr(divide='ignore')
    # channel = np.log(samples)
    channel = np.subtract(channel, zs_LL)
    channel = np.multiply(channel, conv_01_99)
    if invert:
        channel = np.subtract(zg_UL, channel)
    else:
        channel = np.add(zg_LL, channel)

    cmap = loadColormap(palette)
    image = Image.fromarray(np.uint8(cmap(channel)*255))
    
    return image

def loadColormap(palette):
    """ load a Edgetech color map into a matplotlib colormap so  it can be applied to the entire image    """
    df = pd.read_csv(palette, delim_whitespace=True)
    
    # saved_column = df.column_name #you can also use df['column_name']

    cmap = mpl.colors.ListedColormap(['r', 'g', 'b', 'c'])
    
    return cmap
###################################
# zg_LL = lower limit of grey scale
# zg_UL = upper limit of grey scale
# zs_LL = lower limit of samples range
# zs_UL = upper limit of sample range
def samplesToGrayImageLogarithmic(samples, invert, clip):
    zg_LL = 0 # min and max grey scales
    zg_UL = 255
    zs_LL = 0 
    zs_UL = 0
    conv_01_99 = 1
    # channelMin = 0
    # channelMax = 0
    #create numpy arrays so we can compute stats
    channel = np.array(samples)   

    # compute the clips
    if clip > 0:
        channelMin, channelMax = findMinMaxClipValues(channel, clip)
    else:
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
    
    # this scales from the range of image values to the range of output grey levels
    if (zs_UL - zs_LL) is not 0:
        conv_01_99 = ( zg_UL - zg_LL ) / ( zs_UL - zs_LL )
   
    #we can expect some divide by zero errors, so suppress 
    np.seterr(divide='ignore')
    channel = np.log(samples)
    channel = np.subtract(channel, zs_LL)
    channel = np.multiply(channel, conv_01_99)
    if invert:
        channel = np.subtract(zg_UL, channel)
    else:
        channel = np.add(zg_LL, channel)
    # ch = channel.astype('uint8')
    image = Image.fromarray(channel).convert('L')
    
    return image

def angularResponseToPing(ARC, slantRange, numSamples, altitude):
    """compute a pseudo ping of the correct number of samples based on the angular response curve.  This is then easy to apply to the real data""" 
    pseudoPing = [1] * numSamples
    samplesPerMetre = numSamples / slantRange
    altitudeInSamples = math.ceil(altitude * samplesPerMetre)
    if math.isnan(max(ARC)):
        return pseudoPing # nothing to do so quit
    angles = np.arange(0, 90, 1)
    ARC = np.flipud(ARC) #flip the ARC as the numpy interp only works on ascending numbers
    for i in range (altitudeInSamples, numSamples): # skip zero element as it will give div 0 error
        # calculate the angle, then look up the value
        a = math.degrees(math.acos(altitudeInSamples / i))
        pseudoPing[i] = np.interp(a, angles, ARC)
    return pseudoPing
    
def pingToAngularResponse(data, slantRange, numSamples, altitude):
    """Compute an instantaneous angular response using a single ping of information"""
    angularResponse = [1] * 90
    samplesPerMetre = numSamples / slantRange
    altitudeInSamples = altitude * samplesPerMetre
    
    # we do not need to process all records, just subsample on a per-degree basis
    for i in range (0, 90, 1):
        #cos(theta) = A / H
        H = int (altitudeInSamples / math.cos(math.radians(i)))
        if H < len(data): 
            angularResponse[i] = data[H] 
            # print ("angle %d sample %d " % (i, H))    
    # take off angle = arccos (altitude/ samples  
    # angle = math.acos(altitudeInSamples / sample)
    
    return angularResponse     
    
def computeBC(filename, channelA, channelB, samplesPortSum, samplesPortCount, samplesStbdSum, samplesStbdCount, segmentInterval):
    """compute the beam correction table using the sonar ping data for the file. maxSamples lets us know how wide to make the beam correction table.  This is computed by a pre-iteration stage"""
    decimation = 1

    #   open the XTF file for reading 
    print ("Opening file:", filename)
    r = pyXTF.XTFReader(filename)
    while r.moreData():
        ping = r.readPacket()
        segment = altitudeToSegment(ping.SensorPrimaryAltitude, segmentInterval)
        
        # convert the port channel into an array of samples on a per angle basis 
        channelData = np.array(ping.pingChannel[channelA].data[::decimation])
        channelData = np.multiply(channelData, math.pow(2, -ping.pingChannel[channelA].Weight))
        angularResponse = pingToAngularResponse(channelData, ping.pingChannel[channelA].SlantRange, ping.pingChannel[channelA].NumSamples, ping.SensorPrimaryAltitude)
        samplesPortSum[segment] = np.add(samplesPortSum[segment], angularResponse)            
        samplesPortCount[segment] = np.add(samplesPortCount[segment],1)            

        # convert the starboard channel into an array of samples on a per angle basis 
        channelData = np.array(ping.pingChannel[channelB].data[::decimation])
        channelData = np.multiply(channelData, math.pow(2, -ping.pingChannel[channelB].Weight))
        angularResponse = pingToAngularResponse(channelData, ping.pingChannel[channelB].SlantRange, ping.pingChannel[channelB].NumSamples, ping.SensorPrimaryAltitude)
        samplesStbdSum[segment] = np.add(samplesStbdSum[segment], angularResponse)            
        samplesStbdCount[segment] = np.add(samplesStbdCount[segment],1)            
        
        if ping.PingNumber % 1000 == 0:
            print ("Ping: %f" % (ping.PingNumber))                  
    
    print("Complete reading XTF file :-)")
    return (samplesPortSum, samplesPortCount, samplesStbdSum, samplesStbdCount)
    
if __name__ == "__main__":
    main()
