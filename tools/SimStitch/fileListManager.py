import shutil
import tempfile
import os
import optparse
import traceback
import sys
import csv

from galaxy.jobs import JobDestination

import xml.etree.ElementTree as etree
import time, datetime, sys
from xml.dom import minidom

def WriteXML(pars, fnCSV, fnOut):

    root = etree.Element('fileList')
    TimeStamp = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

    root.set('version', "1.0") #XML
    root.set('TimeStamp', TimeStamp) #XML

    ############################################

    setIdentifier = etree.SubElement(root, 'setIdentifier')
    setIdentifier.set("type", "text")

    rootDirectory = etree.SubElement(root, 'rootDirectory')
    rootDirectory.set("type", "text")
    rootDirectory.text=pars[0].replace("__at__", "@")

    rawDirectory = etree.SubElement(root, 'rawDirectory')
    rawDirectory.set("type", "text")

    spec = etree.SubElement(root, 'spec')
    spec.set("type", "struct")

    ############################################

    analysisSplit = etree.SubElement(root, "analysisSplit")
    analysisSplit.set("type", "numeric")
    analysisSplit.text = "0"

    datDirectory = etree.SubElement(root, "datDirectory")
    datDirectory.set("type", "text")


    avgdTransDirectory = etree.SubElement(root, "avgdTransDirectory")
    avgdTransDirectory.set("type", "text")


    overlappingSpecDirectory = etree.SubElement(root, "overlappingSpecDirectory")
    overlappingSpecDirectory.set("type", "text")


    stitchedSpecDirectory = etree.SubElement(root, "stitchedSpecDirectory")
    stitchedSpecDirectory.set("type", "text")


    peakListsDirectory = etree.SubElement(root, "peakListsDirectory")
    peakListsDirectory.set("type", "text")
 

    filteredPeaksDirectory = etree.SubElement(root, "filteredPeaksDirectory")
    filteredPeaksDirectory.set("type", "text")


    blankFlaggedDirectory = etree.SubElement(root, "blankFlaggedDirectory")
    blankFlaggedDirectory.set("type", "text")


    combinedPeaksDirectory = etree.SubElement(root, "combinedPeaksDirectory")
    combinedPeaksDirectory.set("type", "text")

    internalCalFile = etree.SubElement(root, "internalCalFile")
    internalCalFile.set("type", "text")

    numReps = etree.SubElement(root, "numReps")
    numReps.set("type", "numeric")
    numReps.text = str(pars[1])

    notes = etree.SubElement(root, "notes")
    notes.set("type", "text")

    spectraCount = etree.SubElement(root, "spectraCount")
    spectraCount.set("type", "numeric")

    DSO_Directory = etree.SubElement(root, "DSO_Directory")
    DSO_Directory.set("type", "text")

    MultivarModel_Directory = etree.SubElement(root, "MultivarModel_Directory")
    MultivarModel_Directory.set("type", "text")

    MetID_Directory = etree.SubElement(root, "MetID_Directory")
    MetID_Directory.set("type", "text")

    specHeaders = etree.SubElement(root, "specHeaders")
    specHeaders.set("type", "cell")

    specData = etree.SubElement(root, "specData")
    specData.set("type", "cell")

    ##############################################
    with open(fnCSV.replace("__at__", "@"),'rb') as csvfile:

        csvreader = csv.reader(csvfile, delimiter=',',quotechar='"')
    
        count = 0
        
        for row in csvreader:
            if any(row):
                count +=1
                instance = etree.SubElement(spec, 'instance')
                instance.set("index", str(count))
                
                ID = etree.SubElement(instance, 'ID')
                ID.set("type", "text")

                rawFile = etree.SubElement(instance, 'rawFile')
                rawFile.set("type", "text")

                if ".raw" in row[0] or ".RAW" in row[0]:
                    rawFile.text = row[0]
		    ID.text = row[0].rstrip(".raw").rstrip(".RAW")
		    
                else:
                    rawFile.text = "%s.RAW" % (row[0])
		    ID.text = row[0]
                
                sampleID = etree.SubElement(instance, 'sampleID')
                sampleID.set("type", "text")
                sampleID.text = row[1]

    count_xml = etree.SubElement(root, "count")
    count_xml.set("type", "numeric")
    count_xml.text = str(count)

    ##############################################

    out = open(fnOut, "w")

    rough_string = etree.tostring(root, 'utf-8')

    reparsed = minidom.parseString(rough_string)

    #print reparsed.toprettyxml(indent="  ")

    out.write(reparsed.toprettyxml(indent="  "))

    out.close()


def main():
   
    parser = optparse.OptionParser()
    parser.add_option('--outfile',dest='outfile',default='results.xml')
    parser.add_option('--incsv',dest='csvFileList',default=None)
    parser.add_option('--dir', dest='rootDir',default=None )
    parser.add_option('--numreps', dest='numreps', default=3)
	
    (options, args) = parser.parse_args()


    #check for '/' at end of rootDir. It helps.
    if options.rootDir[-1]!='/':
	options.rootDir = options.rootDir+'/'

    #try:
    WriteXML([options.rootDir, options.numreps], options.csvFileList, options.outfile)

    #except:
	#error_file = open(options.outfile, 'w')
	#error_file.write(command_outer)
	#error_file.close()
	#return

    



 

if __name__ == "__main__":
    main()

