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

def WriteXML(pars, fnCSV, fnOut):

    root = etree.Element('fileList')
    TimeStamp = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

    root.set('version', "1.0") #XML
    root.set('TimeStamp', TimeStamp) #XML

    rootDirectory = etree.SubElement(root, 'RootDirectory')
    rootDirectory.set("type", "text")
    rootDirectory.text=pars[0].replace("__at__", "@")

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

                if ".raw" in row[0] or ".RAW" in row[0]:
                    ID.text = row[0].rstrip(".raw").rstrip(".RAW")
                elif ".d" in row[0] or ".D" in row[0]:
                    ID.text = row[0].rstrip(".d").rstrip(".D")  
                elif ".mzML" in row[0] or ".mzml" in row[0] or ".MZML" in row[0]:
                    ID.text = row[0].rstrip(".mzml").rstrip(".mzML").rstrip(".MZML") 
                dataFile.text = row[0]

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
    
    args = parser.parse_args()
    
    if args.dir[-1]!='/':
        args.dir = args.dir+'/'
    
    WriteXML([args.dir, args.numreps], args.incsv, args.outfile)

if __name__ == "__main__":
    main()

