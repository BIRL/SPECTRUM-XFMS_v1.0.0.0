function [Data_file_sort,File_AminoAcidHeader] =  Blocks(ResultsAfterMatch,Oxidize_mz,Unoxidize_mz)
% this funcstion calculates the blocks of data of unoxidized and ozidized mz and place them
% in distinct columns such that the columns from start has Peak, retntention time and area of unxidixzed mz followed by oxidized
%Get the start indices
FileAfterMatch=ResultsAfterMatch;
File_AminoAcidHeader=unique(ResultsAfterMatch(:,1))
File_AminoAcidHeader=File_AminoAcidHeader(2:length(File_AminoAcidHeader))
FileAfterMatch_double_array=str2double(FileAfterMatch);
% to find the index of the file that have m/z values
for Table_index=1:length(File_AminoAcidHeader)
    ind=find(FileAfterMatch(:,1)== string(File_AminoAcidHeader(Table_index)))
    %ind=find(contains(FileAfterMatch(:,1),string(File_AminoAcidHeader(Table_index))))

    %ind=find(contains(file1(:,2),string('TC25')))

    size(ind,1)
    RowCounter = ind( size(ind,1) )+1;
    % extracting those rows that have 0 value
    while ~all(isnan(FileAfterMatch_double_array(RowCounter,1:20))== 1)
        RowCounter = RowCounter + 1;
        endIdx = RowCounter - 1;
    end

    % find the index of unoxidized m/z
    Unoxi_index=find(contains  (ResultsAfterMatch  ( :  ,2),string(Unoxidize_mz)))
    % to find the index of oxidized m/z
    Oxi_index=find(contains  (ResultsAfterMatch (: ,2),string(Oxidize_mz)))
    % defining the range of blocks
    range_lower_limit=ind(1);
    range_upper_limit=ind(length(ind));

    Input_vector=Unoxi_index
    Indices_in_input_vector_within_the_range= intersect(find(Input_vector>=range_lower_limit),find(Input_vector<=range_upper_limit))
    Unoxi_index=Input_vector(Indices_in_input_vector_within_the_range)

    Input_vector2=Oxi_index
    Indices_in_input_vector_within_the_range2= intersect(find(Input_vector2>=range_lower_limit),find(Input_vector2<=range_upper_limit))
    Oxi_index=Input_vector2(Indices_in_input_vector_within_the_range2)
    if isempty(Oxi_index)
        % Extracting the Blocks of the oxidized and un oxidized data
        
        [AAstartId_U,AAendId_U]=UnoxidizeMZCase(Unoxi_index, ind, Oxi_index,FileAfterMatch,Oxidize_mz)

        % rearranging the data such that in the first four columns we have
        % unoxidized data and in the next four columns we have unoxidize data
        Position=1
        while  Position<= length(AAendId_U)
            Data_file_sort(AAstartId_U(Position):AAendId_U(Position),1:4)= ResultsAfterMatch(AAstartId_U(Position):AAendId_U(Position),2:5)
            Position=Position+1;
        end
     else
        % Extracting the Blocks of the oxidized and un oxidized data
        [AAstartId_U,AAendId_U]=UnoxidizeMZ(Unoxi_index, ind, Oxi_index,FileAfterMatch,Oxidize_mz)
        [AAstartId_oxi,AAendId_oxi]=OxidizeMZ(Unoxi_index, ind, Oxi_index,FileAfterMatch,Unoxidize_mz,endIdx)
        % rearranging the data such that in the first four columns we have
        % unoxidized data and in the next four columns we have unoxidize data
        Position=1
        while  Position<= length(AAendId_U)
            Data_file_sort(AAstartId_U(Position):AAendId_U(Position),1:4)= ResultsAfterMatch(AAstartId_U(Position):AAendId_U(Position),2:5)
            Position=Position+1;
        end
            Position=1
            while Position<= length(AAendId_oxi)
                Data_file_sort(AAstartId_oxi(Position):AAendId_oxi(Position),5:19)= ResultsAfterMatch(AAstartId_oxi(Position):AAendId_oxi(Position),2:16)
                Position=Position+1;
            end


    end


end
















