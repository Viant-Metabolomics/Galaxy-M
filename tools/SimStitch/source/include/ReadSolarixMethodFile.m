function out = ReadSolarixMethodFile(xmlfile)

try
	DOMnode = xmlread(xmlfile);
catch
	errordlg('Failed to read XML file.','Import error');
end	

if ~DOMnode.hasChildNodes()
    errordlg('XML tree empty!','Import error');
    return
end

RootNode = DOMnode.getFirstChild;
%RootName = char(RootNode.getNodeName);% capture name of the node

out = [];
if (RootNode.hasChildNodes)
    Child  = RootNode.getChildNodes; % create array of children nodes
    nChild = Child.getLength;
    
    % Read child nodes
    f = [];
    for iChild = 1:nChild        % read in each child
       [cname cLeaf] = NodeName(Child.item(iChild-1), 1);
        if (~isfield(f,cname)),
            f.(cname)=0;           % initialize first time I see this name
        end
        f.(cname) = f.(cname)+1; % add to the counter
        if strcmp(cname,'paramlist')
            ChildNode = Child.item(iChild-1);
            Childsubnodes = ChildNode.getChildNodes;
            nChildsubnodes = Childsubnodes.getLength;
            for iChildsub = 1:nChildsubnodes
                ChildsubNode = Childsubnodes.item(iChildsub-1);
                [csubname, ~] = NodeName(ChildsubNode, 1);
                %[csubname cLeaf] = NodeName(ChildsubNode, 1);
                if ChildsubNode.hasAttributes
                    if strcmp(csubname,'param')
                        %Eerste iChildsub geen content?!
                        Attr = ChildsubNode.getAttributes;
                        Attr_name = char(Attr.item(0).getValue);
                        if f.paramlist == 1
                            counter = 0;
                            if strcmp(Attr_name,'ML1')
                                % Save
                                counter = counter + 1;
                                s = node2struct(ChildsubNode);
                                out.ML1 = str2num(s.children.children.data);
                            elseif  strcmp(Attr_name,'ML2')
                                % Save
                                counter = counter + 1;
                                s = node2struct(ChildsubNode);
                                out.ML2 = str2num(s.children.children.data);
                            elseif  strcmp(Attr_name,'ML3')
                                % Save
                                counter = counter + 1;
                                s = node2struct(ChildsubNode);
                                out.ML3 = str2num(s.children.children.data);
                            elseif strcmp(Attr_name, 'TD_Broadband');
                                counter = counter + 1;
                                s = node2struct(ChildsubNode);
                                out.TD_Broadband = str2num(s.children.children.data);
                            elseif strcmp(Attr_name, 'FR_low');
                                s = node2struct(ChildsubNode);
                                out.FR_low = str2num(s.children.children.data);
                                counter = counter + 1;
                            elseif strcmp(Attr_name, 'SW_h_Broadband');
                                counter = counter + 1;
                                s = node2struct(ChildsubNode);
                                out.SW_h_Broadband = str2num(s.children.children.data);
                            elseif counter == 6
                                break;
                            else
                                continue;
                            end
                        else
                            counter = [];
                            % Get start time and stop time
                            attr = ChildNode.getAttributes;
                            out.segment(f.(cname)-1,1).startTimeMinutes= char(attr.item(0).getValue);
                            out.segment(f.(cname)-1,1).stopTimeMinutes= char(attr.item(1).getValue);
                            if strcmp(Attr_name,'Q1Mass$0')
                                counter = counter + 1;
                                s = node2struct(ChildsubNode);
                                out.segment(f.(cname)-1,1).Q1Mass0=str2num(s.children.children.data);
                            elseif strcmp(Attr_name,'Q1Res$0')
                                counter = counter + 1;
                                s = node2struct(ChildsubNode);
                                out.segment(f.(cname)-1,1).Q1Res0=str2num(s.children.children.data);
                            elseif counter == 2
                                break;
                            end
                        end
                    end
                end
            end
        end
    end
end

