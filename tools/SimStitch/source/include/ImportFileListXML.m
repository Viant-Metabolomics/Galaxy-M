function fileListStruct = ImportFileListXML(file_location)
% 
% function fileListStruct = import_filelist_xml(file_location)
%
% function to allow SimStitch pipeline to read/import xml representations of individual fileList entries 
%
% specifically aimed at the Galaxy Project implementation of SimStitch etc.
% inputs:
%   file_location:  the location of a fileList.xml file as produced by
%                   FileListManagerGUI
%
% outputs:
%   fileListStruct:a matlab struct conforming to SimStitch's fileList
%                   structure
%
% R.L.Davidson 
% 29/05/2014
%
%  J.Engel (11/01/2016): updated code to match output of fileListManager.py
%  script (version from 07/01/2016).

try
    tree = xmlread(file_location);
catch
    errordlg('Failed to read XML file.','Import error');
end

if ~tree.hasChildNodes()
    errordlg('XML tree empty!','Import error');
    return
end

fileListStruct = [];
fileListStruct.Samples= struct('ID',{},'dataFile',{},'sampleID',{});
%fileListStruct.ATDir = '';
fileListStruct.SamplesHeaders = {};
fileListStruct.nReplicates = 0;
fileListStruct.nDataFiles = 0;

%start to interrogate the xml input
if strcmp('fileList',tree.getChildNodes.item(0).getNodeName)
    fileList = tree.getChildNodes.item(0);
else
    errordlg('XML Document root is not fileList - stopping')
    return
end

% if ~isempty(fileList.getElementsByTagName('Instrument').item(0).getChildNodes.item(0))    
%     fileListStruct.Instrument = char(fileList.getElementsByTagName('Instrument').item(0).getChildNodes.item(0).getData);
% end

if ~isempty(fileList.getElementsByTagName('RootDirectory').item(0).getChildNodes.item(0))    
    fileListStruct.RootDirectory = char(fileList.getElementsByTagName('RootDirectory').item(0).getChildNodes.item(0).getData);
end

if ~isempty(fileList.getElementsByTagName('nReplicates').item(0).getChildNodes.item(0))
    fileListStruct.nReplicates = str2num(fileList.getElementsByTagName('nReplicates').item(0).getChildNodes.item(0).getData);
end

if ~isempty(fileList.getElementsByTagName('nDataFiles').item(0).getChildNodes.item(0))
    fileListStruct.nDataFiles = str2num(fileList.getElementsByTagName('nDataFiles').item(0).getChildNodes.item(0).getData);
end


spec_node = fileList.getElementsByTagName('Samples').item(0);
instances = spec_node.getElementsByTagName('instance');

for j = 1:instances.getLength;
    fileListStruct.Samples(j).dataFile = char(instances.item(j-1).getElementsByTagName('dataFile').item(0).getChildNodes.item(0).getData);
    fileListStruct.Samples(j).ID = char(instances.item(j-1).getElementsByTagName('ID').item(0).getChildNodes.item(0).getData);
    if ~isempty(instances.item(j-1).getElementsByTagName('sampleID').item(0).getChildNodes.item(0))
        fileListStruct.Samples(j).sampleID = char(instances.item(j-1).getElementsByTagName('sampleID').item(0).getChildNodes.item(0).getData);
    end
    if ~isempty(instances.item(j-1).getElementsByTagName('batchID').item(0).getChildNodes.item(0))
        fileListStruct.Samples(j).batchID = char(instances.item(j-1).getElementsByTagName('batchID').item(0).getChildNodes.item(0).getData);
    end
    if ~isempty(instances.item(j-1).getElementsByTagName('orderID').item(0).getChildNodes.item(0))
        fileListStruct.Samples(j).orderID = char(instances.item(j-1).getElementsByTagName('orderID').item(0).getChildNodes.item(0).getData);
    end
end

end

