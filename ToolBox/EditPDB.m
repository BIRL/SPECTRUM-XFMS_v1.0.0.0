%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM-XFMS: A MATLAB Toolbox to Analyze X-ray Footprinting Mass Spectrometry Data %
%                                    Version 1.0.0                                     %
%        Copyright (c) Biomedical Informatics & Engineering Research Laboratory,       %
%           Lahore University of Management Sciences Lahore (LUMS), Pakistan.          %
%                           (http://biolabs.lums.edu.pk/BIRL)                          %
%                                (safee.ullah@gmail.com)                               %
%                           Last Modified on: 17-August-2022                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EditPDB(FinalTableOutput,PDBFile)
%This function reads the pdb file and edit it by replacing the temp factor
%values by log(PF) in Chai n A of PDB file 
% pdb File
Pfile = PDBFile;
% final PF_ SASA table 
Tfile=table2cell(FinalTableOutput) ;
Tfile=string(Tfile);

%%%%%% remove the temp values 
for index= 1: length(Pfile.Model.Atom)
    Pfile.Model.Atom(index).tempFactor= '';
end
for index= 1: length(Pfile.Model.HeterogenAtom)
    Pfile.Model.HeterogenAtom(index).tempFactor= '';
end
%%%%% find the position of alpha carbon
if(PDBFile.Model.Atom(1).chainID == 'A')
    pos=find(strcmp({Pfile.Model.Atom(:).AtomName},'CA')& strcmp({Pfile.Model.Atom(:).chainID},'A'));
    pos2=find(strcmp({Pfile.Model.HeterogenAtom(:).AtomName},'CA')& strcmp({Pfile.Model.HeterogenAtom(:).chainID},'A'));
elseif (PDBFile.Model.Atom(1).chainID == 'B')
        pos=find(strcmp({Pfile.Model.Atom(:).AtomName},'CA')& strcmp({Pfile.Model.Atom(:).chainID},'B'));
        pos2=find(strcmp({Pfile.Model.HeterogenAtom(:).AtomName},'CA')& strcmp({Pfile.Model.HeterogenAtom(:).chainID},'B'));
end

% loop to replace the temp values by the PF value
for ind = 1: length(pos)
    j_ind=pos(ind);
ind2= Pfile.Model.Atom(j_ind).resSeq; 
Pfile.Model.Atom(j_ind).tempFactor = Tfile(ind2-1,4);
end
% loop to repalce missing values with blanks 
for ind= 1: length(Pfile.Model.Atom)
if ismissing( Pfile.Model.Atom(ind).tempFactor)
    Pfile.Model.Atom(ind).tempFactor = '';
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in case of heterogenous atoms 
for ind = 1: length(pos2)
    j_ind=pos2(ind);
ind2= Pfile.Model.HeterogenAtom(j_ind).resSeq; 
Pfile.Model.HeterogenAtom(j_ind).tempFactor = Tfile(ind2-1,4);
end

for ind= 1: length(Pfile.Model.HeterogenAtom)
if ismissing( Pfile.Model.HeterogenAtom(ind).tempFactor)
    Pfile.Model.HeterogenAtom(ind).tempFactor = '';
end
end

cd ('Result')
pdbwrite('Modifiedchey.pdb',Pfile)
