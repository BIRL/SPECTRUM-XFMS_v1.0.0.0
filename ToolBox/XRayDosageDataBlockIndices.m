% This function generate blocks of data according to the dosage Values
function [row_ind, row_end] =  XRayDosageDataBlockIndices(file)
%Get the start indices
row_ind = find(cellfun('length',regexp(string(file),'.x')) == 1);
row_end = [];
%cahnge the string file to double
file = str2double(file);
FileTotalRows=size(file)
FileTotalRows=FileTotalRows(1)
%Get the end indices
for i = 1: size(row_ind,1)
    %start index is equal to the row value at position i
    startIdx = row_ind(i);
    if i == size(row_ind,1)
        %RowCounter=length(file)
       % RowCounter = row_ind( size(row_ind,1) )+1;
        % extracting those rows that have 0 value
        %while ~all(isnan(file(RowCounter,:))== 1)
           % RowCounter = RowCounter + 1;
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
