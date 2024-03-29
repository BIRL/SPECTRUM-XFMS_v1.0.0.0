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
function [row_ind, row_end] =  XRayDosageDataBlockIndices(file)
%Get the start indices
row_ind = find(cellfun('length',regexp(string(file),'.x')) == 1);
row_end = [];
%change the string file to double
file = str2double(file);
FileTotalRows=size(file);
FileTotalRows=FileTotalRows(1);
%Get the end indices
for i = 1: size(row_ind,1)
    %start index is equal to the row value at position i
    startIdx = row_ind(i);
    if i == size(row_ind,1)
                  endIdx = FileTotalRows;
        %end
    else
        % end index is defined as the value at i+1 in row_ind
        endIdx = row_ind(i+1)-1; 
    end
    RowData = file(endIdx-1,:);
    while isnan(RowData(:,:))
        % assigning end id as the endIdx - 1
        endIdx = endIdx - 1;  
        RowData = file(endIdx,:);
    end
    row_end = [row_end; endIdx];
end
