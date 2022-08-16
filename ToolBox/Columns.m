%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM-XFMS: A MATLAB Toolbox to Analyze X-ray Footprinting Mass Spectrometry Data %
%                                    Version 1.0.0                                     %
%        Copyright (c) Biomedical Informatics & Engineering Research Laboratory,       %
%           Lahore University of Management Sciences Lahore (LUMS), Pakistan.          %
%                           (http://biolabs.lums.edu.pk/BIRL)                          %
%                                (safee.ullah@gmail.com)                               %
%                            Last Modified on: 26-July-2022                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Index_var, colname] =  Columns(mock)
%Fine cells that are empty
Index_var=find( cellfun( @isempty, mock(2,:) ) ==0);
%Storage variable for non-empty column names
colname=string({ });
%For all columns in the first row that are not empty
for i=1: length(Index_var)
    %Read data from file
    ColumnInfo= mock(2,Index_var(i));
    %Store column name
    colname(1,Index_var(i))=ColumnInfo;
end
end