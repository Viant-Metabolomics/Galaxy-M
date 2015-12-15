% The code below was taken from the xml_read function by Jaroslaw Tuszynski
% Copyright (c) 2007
function str = str2varName(str, KeepNS)
% convert a sting to a valid matlab variable name
if(KeepNS)
  str = regexprep(str,':','_COLON_', 'once', 'ignorecase');
else
  k = strfind(str,':');
  if (~isempty(k))
    str = str(k+1:end);
  end
end
str = regexprep(str,'-','_DASH_'  ,'once', 'ignorecase');
if (~isvarname(str)) && (~iskeyword(str))
  str = genvarname(str);
end

function [Name LeafNode] = NodeName(node, KeepNS)
% get node name and make sure it is a valid variable name in Matlab.
% also get node type:
%   LeafNode=0 - normal element node,
%   LeafNode=1 - text node
%   LeafNode=2 - supported non-text leaf node,
%   LeafNode=3 - supported processing instructions leaf node,
%   LeafNode=-1 - unsupported non-text leaf node
switch (node.getNodeType)
  case node.ELEMENT_NODE
    Name = char(node.getNodeName);% capture name of the node
    Name = str2varName(Name, KeepNS);     % if Name is not a good variable name - fix it
    LeafNode = 0;
  case node.TEXT_NODE
    Name = 'CONTENT';
    LeafNode = 1;
  case node.COMMENT_NODE
    Name = 'COMMENT';
    LeafNode = 2;
  case node.CDATA_SECTION_NODE
    Name = 'CDATA_SECTION';
    LeafNode = 2;
  case node.DOCUMENT_TYPE_NODE
    Name = 'DOCUMENT_TYPE';
    LeafNode = 2;
  case node.PROCESSING_INSTRUCTION_NODE
    Name = 'PROCESSING_INSTRUCTION';
    LeafNode = 3;
  otherwise
    NodeType = {'ELEMENT','ATTRIBUTE','TEXT','CDATA_SECTION', ...
      'ENTITY_REFERENCE', 'ENTITY', 'PROCESSING_INSTRUCTION', 'COMMENT',...
      'DOCUMENT', 'DOCUMENT_TYPE', 'DOCUMENT_FRAGMENT', 'NOTATION'};
    Name = char(node.getNodeName);% capture name of the node
    warning('xml_io_tools:read:unkNode', ...
      'Unknown node type encountered: %s_NODE (%s)', NodeType{node.getNodeType}, Name);
    LeafNode = -1;
end

