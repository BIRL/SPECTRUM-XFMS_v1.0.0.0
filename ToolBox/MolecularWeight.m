%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM-XFMS: A MATLAB Toolbox to Analyze X-ray Footprinting Mass Spectrometry Data %
%                                    Version 1.0.0                                     %
%        Copyright (c) Biomedical Informatics & Engineering Research Laboratory,       %
%           Lahore University of Management Sciences Lahore (LUMS), Pakistan.          %
%                           (http://biolabs.lums.edu.pk/BIRL)                          %
%                                (safee.ullah@gmail.com)                               %
%                            Last Modified on: 26-July-2022                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Oxidize_mz, Unoxidize_mz, Theoretical_peptide_weight] = MolecularWeight( Sequence )
% This function calculates the monoisotopic mass of peptide and compute the
% oxidized AND UNOXIDIZED M/Z
% Sequence is the sequnenc of peptide under consideration
%Sequence='FLVVDDFSTMR'
sequence = char(Sequence);
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];
Average_amino_acid_masses = [71.03779, 156.1857, 114.1026, 115.0874, 103.1429, 129.114, 128.1292, 57.0513, 137.1393, 113.1576, 113.1576,...
    128.1723, 131.1961, 147.1739, 97.1152, 87.0773, 101.1039, 186.2099, 163.1733, 99.1311];
Monoisotopic_amino_acid_masses = [71.03711, 156.10111, 114.04293, 115.02694, 103.00919, 129.04259, 128.05858, 57.02146, 137.05891, 113.08406, 113.08406,...
    128.09496, 131.04049, 147.06841, 97.05276, 87.03203, 101.04768, 186.07931, 163.06333, 99.06841];
[~, sizeofsequence] = size(sequence);

Theoretical_peptide_weight = 0;

% along the complete length of sequence
for alongsequencesize = 1:sizeofsequence
    for alongaminoacids = 1:20

        if amino_acids(alongaminoacids)== sequence(alongsequencesize)
           Theoretical_peptide_weight = Theoretical_peptide_weight + Average_amino_acid_masses(alongaminoacids);
        end
    end
end
Theoretical_peptide_weight = Theoretical_peptide_weight + 18.01524;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  unoxidize m/z  WITH CHARGE STATE 2 nd 3
IndexTable=1;
Unoxidize= (Theoretical_peptide_weight + 2)/2;
Unoxidize_mz1(IndexTable)= string(floor(Unoxidize));
Unoxidize= (Theoretical_peptide_weight + 3)/3;
IndexTable=IndexTable+1;
Unoxidize_mz1(IndexTable) = string(floor(Unoxidize));
Unoxidize= (Theoretical_peptide_weight + 4)/4;
IndexTable=IndexTable+1;
Unoxidize_mz1(IndexTable) = string(floor(Unoxidize));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  oxidize m/z 
%%%%%%%%%%% OXIDATION NUMBER = 2,3,4 and charge state = 2 3
IndexTable=1;
oxidize= (Theoretical_peptide_weight + 2+16)/2;
Oxidize_mz1(IndexTable)= string(floor(oxidize));
oxidize= (Theoretical_peptide_weight + 3+16)/3;
IndexTable=IndexTable+1;
Oxidize_mz1(IndexTable) = string(floor(oxidize));
oxidize= (Theoretical_peptide_weight + 4+16)/4;
IndexTable=IndexTable+1;
Oxidize_mz1(IndexTable) = string(floor(oxidize));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oxidize= (Theoretical_peptide_weight + 2+32)/2;
IndexTable=IndexTable+1;
Oxidize_mz1(IndexTable) = string(floor(oxidize));
oxidize= (Theoretical_peptide_weight + 3+32)/3;
IndexTable=IndexTable+1;
Oxidize_mz1(IndexTable) = string(floor(oxidize));
oxidize= (Theoretical_peptide_weight + 4+32)/4;
IndexTable=IndexTable+1;
Oxidize_mz1(IndexTable) = string(floor(oxidize));
oxidize= (Theoretical_peptide_weight + 2+48)/2;
IndexTable=IndexTable+1;
Oxidize_mz1(IndexTable)= string(floor(oxidize));
oxidize= (Theoretical_peptide_weight + 3+48)/3;
IndexTable=IndexTable+1;
Oxidize_mz1(IndexTable) = string(floor(oxidize));
oxidize= (Theoretical_peptide_weight + 4+48)/4;
IndexTable=IndexTable+1;
Oxidize_mz1(IndexTable) = string(floor(oxidize));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Tolernce for Molecular  weights
IndexTabforTolernce=1;
tol= 1;

for IndexTable= 1:length(Oxidize_mz1)
    Oxidize_mz(IndexTabforTolernce)= str2double(Oxidize_mz1(IndexTable))+tol;
    IndexTabforTolernce=IndexTabforTolernce+1;
     Oxidize_mz(IndexTabforTolernce)=  str2double(Oxidize_mz1(IndexTable))-tol;
     IndexTabforTolernce=IndexTabforTolernce+1;
     Oxidize_mz(IndexTabforTolernce)= Oxidize_mz1(IndexTable);
     IndexTabforTolernce=IndexTabforTolernce+1;

end
IndexTabforTolernce=1;
tol= 1;
for IndexTable= 1:length(Unoxidize_mz1)
    Unoxidize_mz(IndexTabforTolernce)= str2double(Unoxidize_mz1(IndexTable))+tol;
    IndexTabforTolernce=IndexTabforTolernce+1;
     Unoxidize_mz(IndexTabforTolernce)=  str2double(Unoxidize_mz1(IndexTable))-tol;
     IndexTabforTolernce=IndexTabforTolernce+1;
     Unoxidize_mz(IndexTabforTolernce)= Unoxidize_mz1(IndexTable)
     IndexTabforTolernce=IndexTabforTolernce+1;

end



