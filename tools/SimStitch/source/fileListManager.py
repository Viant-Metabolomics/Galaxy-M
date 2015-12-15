import argparse
import shutil
import tempfile
import os
import traceback
import sys
import csv

import xml.etree.ElementTree as etree
import time, datetime, sys
from xml.dom import minidom

def WriteXML(pars, instrument, fnCSV, fnOut):

    root = etree.Element('fileList')
    TimeStamp = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

    root.set('version', "1.0") #XML
    root.set('TimeStamp', TimeStamp) #XML

    rootDirectory = etree.SubElement(root, 'RootDirectory')
    rootDirectory.set("type", "text")
    rootDirectory.text=pars[0].replace("__at__", "@")

    instrumentType = etree.SubElement(root, 'Instrument')
    instrumentType.set("type", "text")
    instrumentType.text = instrument

    spec = etree.SubElement(root, 'Samples')
    spec.set("type", "struct")

    with open(fnCSV.replace("__at__", "@"),'rb') as csvfile:

        csvreader = csv.reader(csvfile, delimiter=',',quotechar='"')
    
        c = 0
        
        for row in csvreader:
            if any(row):

                c +=1

                instance = etree.SubElement(spec, 'instance')
                instance.set("index", str(c))
                
                ID = etree.SubElement(instance, 'ID')
                ID.set("type", "text")

                dataFile = etree.SubElement(instance, 'dataFile')
                dataFile.set("type", "text")

                if instrument == "ltqft" or instrument == "orbitrap" or instrument == "qexactive":
                    if ".raw" in row[0] or ".RAW" in row[0]:
                        dataFile.text = row[0]
                        ID.text = row[0].rstrip(".raw").rstrip(".RAW")

                    else:
                        dataFile.text = "%s.RAW" % (row[0])
                        ID.text = row[0]

                elif instrument == "solarix":
                    if ".d" in row[0] or ".D" in row[0]:
                        dataFile.text = row[0]
                        ID.text = row[0].rstrip(".d").rstrip(".D")
                        
                    else:
                        dataFile.text = "%s.d" % (row[0])
                        ID.text = row[0]
                
                sampleID = etree.SubElement(instance,'sampleID')
                sampleID.set("type", "text")
                sampleID.text = row[1]

                batchID = etree.SubElement(instance,'batchID')
                batchID.set("type", "text")
                batchID.text = row[2]

                orderID = etree.SubElement(instance,'orderID')
                orderID.set("type","text")
                orderID.text = row[3]

    nDataFiles = etree.SubElement(root, "nDataFiles")
    nDataFiles.set("type", "numeric")
    nDataFiles.text = str(c)

    nReplicates = etree.SubElement(root, "nReplicates")
    nReplicates.set("type", "numeric")
    nReplicates.text = str(pars[1])

    ##############################################

    out = open(fnOut, "w")

    rough_string = etree.tostring(root, 'utf-8')

    reparsed = minidom.parseString(rough_string)

    out.write(reparsed.toprettyxml(indent="  "))

    out.close()


def main():
    
    parser = argparse.ArgumentParser(description='Create filelist')
    parser.add_argument('--incsv', help='CSV List',required=True)
    parser.add_argument('--outfile',help='Output XML file', required=True, default='filelist.xml')
    parser.add_argument('--dir', help='Directory of Data files', default=None)
    parser.add_argument('--numreps', help='Number of technical replicates', default=3, type = int)
    parser.add_argument('--instrument',help='Instrument Name (QExactive; Orbitrap or LTQFT; or Solarix)', required=True)
    
    args = parser.parse_args()
    
    if args.instrument.lower() not in ["qexactive", "orbitrap", "ltqft", "solarix"]:
        print "Instrument name not available"

    if args.dir[-1]!='/':
        args.dir = args.dir+'/'
    
    WriteXML([args.dir, args.numreps], args.instrument.lower(), args.incsv, args.outfile)

if __name__ == "__main__":
    main()

