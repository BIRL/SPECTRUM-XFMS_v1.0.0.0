%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM-XFMS: A MATLAB Toolbox to Analyze X-ray Footprinting Mass Spectrometry Data %
%                                    Version 1.0.0                                     %
%        Copyright (c) Biomedical Informatics & Engineering Research Laboratory,       %
%           Lahore University of Management Sciences Lahore (LUMS), Pakistan.          %
%                           (http://biolabs.lums.edu.pk/BIRL)                          %
%                                (safee.ullah@gmail.com)                               %
%                            Last Modified on: 26-July-2022                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  MassToChargeRatioList(wholeSeq)
AminoAcidSeq=wholeSeq;
Cleave_parts = cleave(AminoAcidSeq,'trypsin');
MASS=cell(100, 10);
for i=1: size(Cleave_parts,1)

sequence = char(Cleave_parts(i));
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
MASS(i,1)=Cleave_parts(i);
MASS(i,2)= {Theoretical_peptide_weight};
%%%%%%%%% calculate m/z
MASS(i,3)= {(Theoretical_peptide_weight + 2)/2};
MASS(i,4)= {(Theoretical_peptide_weight + 2+16)/2};
MASS(i,5)= {(Theoretical_peptide_weight + 2+32)/2};
MASS(i,6)= {(Theoretical_peptide_weight + 2+48)/2};
%%%%%%%%%%%%%%%
MASS(i,7)= {(Theoretical_peptide_weight + 3)/3};
MASS(i,8)= {(Theoretical_peptide_weight + 3+16)/3};
MASS(i,9)= {(Theoretical_peptide_weight + 3+32)/3};
MASS(i,10)= {(Theoretical_peptide_weight + 3+48)/3};
%%%%%%%%%%%%%%
MASS(i,11)= {(Theoretical_peptide_weight + 4)/4};
MASS(i,12)= {(Theoretical_peptide_weight + 4+16)/4};
MASS(i,13)= {(Theoretical_peptide_weight + 4+32)/4};
MASS(i,14)= {(Theoretical_peptide_weight + 4+48)/4};

end
header= {'peptideSeq','PeptideMass','m/z +2','m/z +2+1*','m/z +2+2*','m/z +2+3*','m/z +3','m/z +3+1*','m/z +3+2*','m/z +3+3*','m/z +4','m/z +4+1*','m/z +4+2*','m/z +4+3*'};
baseFileName='MasstoChargeRatio.xls';
OutputDir = uigetdir(pwd,'Select the Output folder' );
fullFileName = fullfile(OutputDir, baseFileName);
xlswrite(fullFileName,header,'Sheet1','A1');
xlswrite(fullFileName,MASS,'Sheet1','A2');

