#!/usr/bin/python
import argparse
import sys, os
from comtypes.client import CreateObject 
from numpy import array, zeros, object as nobj
from scipy.io import savemat 
from ctypes import c_long, c_double
from comtypes.automation import VARIANT, BSTR
from collections import OrderedDict

def GetRawHandle(RawFile):
    hRaw = CreateObject('MSFileReader.XRawfile')
    hRaw.open(RawFile)
    return hRaw

def RetrieveData(hRaw):

    MatPyLib = {
        'Reply':0,
        'DetectorCount':0,
        'ScanList': [],
        'Filters': [],
        'FiltersCount':0,
        'ScanInfo': OrderedDict()
    }

    pnNumControllers = c_long()
    hRaw.GetNumberOfControllers(pnNumControllers)
    MatPyLib['DetectorCount'] = pnNumControllers.value
    if pnNumControllers.value == 0 or pnNumControllers.value > 1:
        print "Warning: Number of Controlers is 0 or larger than 1"

    nControllerType = 0 #MS
    pnNumControllers = c_long()
    hRaw.GetNumberOfControllersOfType(nControllerType, pnNumControllers)
    if pnNumControllers.value > 1:
        print "Warning: Multipe ControllersOfType of MS"

    reply = 0
    hRaw.SetCurrentController(0,1)

    pnFirstSpectrum = c_long()
    pnLastSpectrum = c_long()
    hRaw.GetFirstSpectrumNumber(pnFirstSpectrum)
    hRaw.GetLastSpectrumNumber(pnLastSpectrum)

    MatPyLib['ScanList'] = range(pnFirstSpectrum.value,pnLastSpectrum.value+1)

    for i in MatPyLib['ScanList']:

        ScanNumber = "Scan" + str(int(i))

        nScanNumber = c_long(int(i))
        
        MatPyLib['ScanInfo'][ScanNumber] = {}

        pbstrFilter = BSTR()
        hRaw.GetFilterForScanNum(nScanNumber, pbstrFilter)
        MatPyLib['ScanInfo'][ScanNumber]["Filter"] = pbstrFilter.value

        pnNumPackets = c_long()
        pdStartTime = c_double()
        pdLowMass = c_double()
        pdHighMass = c_double()
        pdTIC  = c_double()
        pdBasePeakMass  = c_double()
        pdBasePeakIntensity  = c_double()
        pnNumChannels = c_long()
        pbUniformTime = c_long()
        pdFrequency = c_double()

        hRaw.GetScanHeaderInfoForScanNum(nScanNumber,
        pnNumPackets,
        pdStartTime, 
        pdLowMass,
        pdHighMass,
        pdTIC,
        pdBasePeakMass,
        pdBasePeakIntensity,
        pnNumChannels,
        pbUniformTime,
        pdFrequency)

        MatPyLib['ScanInfo'][ScanNumber]["TIC"] = float(pdTIC.value)
        MatPyLib['ScanInfo'][ScanNumber]["LowMass"] = float(pdLowMass.value)
        MatPyLib['ScanInfo'][ScanNumber]["HighMass"] = float(pdHighMass.value)

        pvarLabels = VARIANT()
        pvarValues = VARIANT()
        pnArraySize = c_long()
        
        hRaw.GetTrailerExtraForScanNum(nScanNumber, pvarLabels, pvarValues, pnArraySize)

        MatPyLib['ScanInfo'][ScanNumber]["IT"] = float(pvarValues.value[pvarLabels.value.index('Ion Injection Time (ms):')])
        MatPyLib['ScanInfo'][ScanNumber]["A"] = float(pvarValues.value[pvarLabels.value.index('Conversion Parameter A:')])
        MatPyLib['ScanInfo'][ScanNumber]["B"] = float(pvarValues.value[pvarLabels.value.index('Conversion Parameter B:')])

        szFilter = BSTR()
        nIntensityCutoffType = c_long(0)
        nIntensityCutoffValue = c_long(0)
        nMaxNumberOfPeaks = c_long(0)
        bCentroidResult = c_long(0)
        pdCentroidPeakWidth = c_double(0)
        pvarMassList = VARIANT()
        pvarPeakFlags = VARIANT()
        pnArraySize = c_long()
        
        hRaw.GetMassListFromScanNum(nScanNumber,
        szFilter,
        nIntensityCutoffType,
        nIntensityCutoffValue,
        nMaxNumberOfPeaks,
        bCentroidResult,
        pdCentroidPeakWidth,
        pvarMassList,
        pvarPeakFlags, pnArraySize)

        MatPyLib['ScanInfo'][ScanNumber]["scanmz"] = array(pvarMassList.value[0])
        MatPyLib['ScanInfo'][ScanNumber]["scand"] = array(pvarMassList.value[1])
        
    pvarFilterArray = VARIANT()
    pnArraySize = c_long()
    hRaw.GetFilters(pvarFilterArray, pnArraySize)

    MatPyLib['FiltersCount'] = pnArraySize.value
    MatPyLib["Filters"] = zeros((pnArraySize.value,), dtype = nobj)
    MatPyLib["Filters"][:] = pvarFilterArray.value
    return MatPyLib

def main():
    parser = argparse.ArgumentParser(description='This is a python implemantion of MSFilReader by Ralf Weber.')
    parser.add_argument('-i','--input', help='Input RAW file',required=True)
    parser.add_argument('-o','--output',help='Output Mat file', required=True)
    args = parser.parse_args()
    if os.path.isfile(args.input):
        hRaw = GetRawHandle(args.input)
        DataRaw = RetrieveData(hRaw)
        savemat(args.output, {"temp" :DataRaw}, do_compression=True)
    else:
        print "RAW file does not exist"

    current_pid = os.getpid()
    os.system("taskkill /pid %s /f" % current_pid)

if __name__ == "__main__":
    main()


 

