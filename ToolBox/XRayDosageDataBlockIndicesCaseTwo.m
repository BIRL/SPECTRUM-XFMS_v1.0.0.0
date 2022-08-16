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
% This will generate index of blocks of data of File Two in case of same or
% overlapping header.
function [row_ind2, row_end2] =  XRayDosageDataBlockIndicescaseTwo(FileMat2)
%Get the start indices
row_ind2 = find(cellfun('length',regexp(string(FileMat2),'.x')) == 1);
row_end2 = [];
%cahnge the string file to double
FileMat2 = str2double(FileMat2);
FileTotalRows=size(FileMat2);
FileTotalRows=FileTotalRows(1);
%Get the end indices
for i = 1: size(row_ind2,1)
    %start index is equal to the row value at position i
    startIdx = row_ind2(i);
    if i == size(row_ind2,1)
                endIdx = FileTotalRows;
     
    else
        % end index is defined as the value at i+1 in row_ind
        endIdx = row_ind2(i+1)-1; 
    end
    RowData = FileMat2(endIdx-1,:);
    
    row_end2 = [row_end2; endIdx];
end
