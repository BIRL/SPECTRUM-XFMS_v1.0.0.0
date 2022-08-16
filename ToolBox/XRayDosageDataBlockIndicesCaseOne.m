%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM-XFMS: A MATLAB Toolbox to Analyze X-ray Footprinting Mass Spectrometry Data %
%                                    Version 1.0.0                                     %
%        Copyright (c) Biomedical Informatics & Engineering Research Laboratory,       %
%           Lahore University of Management Sciences Lahore (LUMS), Pakistan.          %
%                           (http://biolabs.lums.edu.pk/BIRL)                          %
%                                (safee.ullah@gmail.com)                               %
%                            Last Modified on: 26-July-2022                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generate blocks of data according to the dosage Values
% This function is only used in case if two peptides have same or
% overlapping headers.
% This will generate index of blocks of data of File One in case of same or
% overlapping header.
function [row_ind1, row_end1] =  XRayDosageDataBlockIndicesCaseOne(FileMat)
%Get the start indices
row_ind1 = find(cellfun('length',regexp(string(FileMat),'.x')) == 1);
row_end1 = [];
%cahnge the string file to double
FileMat = str2double(FileMat);
FileTotalRows=size(FileMat);
FileTotalRows=FileTotalRows(1);
%Get the end indices
for i = 1: size(row_ind1,1)
    %start index is equal to the row value at position i
    startIdx = row_ind1(i);
    if i == size(row_ind1,1)
                  endIdx = FileTotalRows;
   
    else
        % end index is defined as the value at i+1 in row_ind
        endIdx = row_ind1(i+1)-1; 
    end
    
    row_end1 = [row_end1; endIdx];
end
