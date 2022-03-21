%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           SPECTRUM-XFMS                          %
%                           Version 1.0.0.0                        %
%        Copyright (c) Biomedical Informatics Research Laboraty,   %
%          Lahore University of Management Sciences Lahore (LUMS), %
%                           Pakistan.                              %
%                (http://biolabs.lums.edu.pk/BIRL)                 %
%                    (safee.ullah@gmail.com)                       %
%                 Last Modified on: 20-March-2022                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MainCalculation()
% This function reads the following Inputs:
% 1. folder containing mass hunter files of each peptide
% 2. Mascot file
% 3. Fasta file of protein
% It matches the mass hunter files with the mascost file on basis of
% Retention time and give us a generalized output for each replicate of a
% peptide that can be used as an input for the 2nd part (dose response)
% The generalized output file contains following information:
% 1: X-ray dosage
% 2: Unxodized area
% 3: Oxdized area
% 4: Peak split  information
% 5: Oxidized residues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step1 : Load the Mascot , masshunter, SASA , Protein fasta File
% Get the directory from user and read the Mascot file
InputDir = uigetdir(pwd,'Select the Input folder' )
MascotDir = uigetfile({'*xls;*.fasta;*.xlsx'},'Select a Mascot File' )
%Mfiles = dir(fullfile(MascotDir, '*.xlsx'))
dir_mascot= char(strcat(InputDir,strcat('\',MascotDir)))

%Reading mascot file
[~,~,mascotfile] = xlsread(dir_mascot);
mascotfile = string(mascotfile);

%%%%%% Reading the FASTA file
%select the fasta file
FASTADir = uigetfile({'*xls;*.fasta;*.xlsx'},'Select a FASTA File' )
dir_fasta= char(strcat(InputDir,strcat('\',FASTADir)))
% Read fasta file
wholeSeq=fastaread(dir_fasta)
% geting the sequence from the fasta file
wholeSeq=wholeSeq.Sequence
%%%%%%%%%%%%%%% Seading SASA file
SASADir = uigetfile({'*xls;*.fasta;*.xlsx'},'Select a SASA File' )
dir_SASA= char(strcat(InputDir,strcat('\',SASADir)))

FileSASA=readtable(dir_SASA);
FileSASA = table2cell(FileSASA)
%%%%%%%%%%%%%%%%%% Reading MassHunter  files
ProjectData = uigetdir(pwd,'Select a folder which contains MassHunter Data files' );

files = dir(ProjectData);

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];

% Extract only those that are directories.
Peptides = files(dirFlags); % A structure with extra info.

AllPeptides = {};
for number=1:size(Peptides,1)
    CurrentPeptide = Peptides(number).name;
    if(strcmp (CurrentPeptide, '.') || isempty(CurrentPeptide) || strcmp (CurrentPeptide, '..'))
        %do nothing
    else
        AllPeptides = [AllPeptides; CurrentPeptide];

        %1. Go into directory of Peptide and find the sub directories of
        %  Relicates
        dir_sub= char(strcat(ProjectData,strcat('\',Peptides(number).name)))
        Replicates= dir(dir_sub)
        % Read the replicates directory one after other
        AllReplicates = {};
        for num=1:size(Replicates,1)
            CurrentReplicate = Replicates(num).name;
            if (strcmp (CurrentReplicate, '.') || isempty(CurrentReplicate) || strcmp (CurrentReplicate, '..'))
                %do nothing
            else
                AllReplicates = [AllReplicates; CurrentReplicate];
                dir_rep= char(strcat(dir_sub,strcat('\',Replicates(num).name)))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                MassHunterDataFiles=dir(fullfile(dir_rep,'\*.xls'))
                % Get a list of all files in the folder with the desired file name pattern.
                MassHunterDataFiles=natsortfiles(MassHunterDataFiles)
                %Tolerance used for matching Peaks
                MatchedTolerance = 0.06;
                size_of_mascot = size(mascotfile);
                num_of_files = length(MassHunterDataFiles);

                for i=1:size(mascotfile,2)
                    %Takes column number of m/z values (887, 895, 592)
                    if mascotfile(1,i)== 'pep_exp_mz'
                        pep_exp_mz_column = i;
                        %Stores the column number of RT
                    elseif mascotfile(1,i)== 'pep_scan_title'
                        pep_scan_title_column = i;
                    elseif mascotfile(1,i)== 'pep_var_mod_pos'
                        pep_var_mod_pos_column = i;
                        % store the peptide sequence
                    elseif mascotfile(1,i)== 'pep_seq'
                        pep_seq_column = i;
                    end
                end
                % defining the RT column and Area column of mass Hunter files
                RT_column = 2;
                Area_column = 3;
                uniqueheader = ["Peak", "RT", "Area", "Match1", "Residue1", "Match2", "Residue2", "Match3", "Residue3", "Match4", "Residue4", "Match5", "Residue5"];

                %Data for Mega Matrix for Handling Result Excel File Data
                MegaMatrix = strings(300,50);
                % this index will help in populating MegaMatrix w.r.t. X-axis
                indexMegaMatrixXaxis = 1;
                % this index will help in populating MegaMatrix w.r.t. Y-axis
                indexMegaMatrixYaxis = 8;
                % Number of .xlsx files present in a folder
                for indexMassHunter=1:num_of_files(1)
                    FlagSetforHeader = 0;

                    subindexStartforArea = 0;
                    subindexforArea = 0;

                    %indexMegaMatrixXaxis = indexMegaMatrixXaxis + 1;
                    % Step 2: Read the MassHunter files
                    num = 0;
                    File = MassHunterDataFiles(indexMassHunter).name;
                    % Extracting the complete name (Location + Name)
                    FullFileName = fullfile(MassHunterDataFiles(indexMassHunter).folder, File);
                    % Read the .xlsx file
                    [num, txt, MassHunterDataFile] = xlsread(FullFileName);

                    sizeofMassHunterDataFile = size(MassHunterDataFile)
                    MassHunterDataFile = string(MassHunterDataFile);
                    MassHunterDataFile= MassHunterDataFile(:,1:3)
                    %Extracting mz value from excel file name
                    mzvalue_indecimal=extractBetween(File,'_','.xls');
                    %Iterations are also for progressbar
                    iteration = indexMassHunter;
                    %Dose value
                    dosevalue=extractBetween(File,1,'_')
                    header = [string(dosevalue),string(dosevalue)+" "+ string(mzvalue_indecimal)];
                    count=0;
                    %Removes the decimal part of mz value %#FUTURECOMPATIBILITY
                    mzvalue = string(floor(str2double(mzvalue_indecimal)))
                    %%Getting Rows where specific mzvalue Start and Ends into the Mascot file e.g. 597 starts at 25 and ends at 31

                    % for start index
                    for index=1:size_of_mascot(1)
                        mz1=string(mascotfile(index,pep_exp_mz_column))
                        %Taking mz value from pep_exp_mz column
                        if extractBetween(mz1,1,'.')==mzvalue
                            count_start=index;
                            break
                        end
                    end
                    for itermascotfiledata=1:size_of_mascot(1)
                        mz1=string(mascotfile(itermascotfiledata,pep_exp_mz_column));
                        %Takes the number of times an mz value is repeated in mascot file
                        if extractBetween(mz1,1,'.')==mzvalue
                            count=count+1;
                        end
                    end
                    % For Compensating last value of mascot file
                    if count ~= 1
                        count_total=count_start+count-1;
                    else
                        count_start;
                        count=1;
                        count_total = count_start;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %step 3: Match the MassHunter Rt with the Mascot RT  with a tolerence of 0.06
                    % Scaning RT values and Area % Scaning from 2nd Row, 1st Row is Header
                    for indexMassHunterFileData=2:sizeofMassHunterDataFile(1)
                        % We want to write just one time Mass Hunter Peak Info e.g. (Peak RT Area Height Type Width FWHM) so its a flag
                        peakchange = 0;
                        indexMegaMatrixYaxis = 8;
                        RT_to_be_matched=MassHunterDataFile(indexMassHunterFileData,RT_column);
                        if ~ismissing(RT_to_be_matched)
                            Area_of_RT=str2double(MassHunterDataFile(indexMassHunterFileData,Area_column));
                            %%Defining tolerance of +/- 0.06
                            upper_limit=str2double(RT_to_be_matched)+MatchedTolerance;
                            lower_limit=str2double(RT_to_be_matched)-MatchedTolerance;
                            % Extract pep_exp_mz from Mascot file e.g. 592 etc.
                            for j=count_start:count_total
                                Sequence=char(mascotfile(j,pep_seq_column))
                                peak_str=extractBetween(mascotfile(j,pep_scan_title_column),'at ',' mins');
                                peak=str2double(peak_str);
                                if lower_limit<peak && peak<upper_limit
                                    RT=str2double(RT_to_be_matched);
                                    % This is for UnOxidized Peak Information      (MAKE SEPRATE FUNCTION FOR FORMATING)
                                    if ismissing(mascotfile(j,pep_var_mod_pos_column))

                                        if FlagSetforHeader == 0
                                            MegaMatrix(indexMegaMatrixXaxis, 1:2) = header;

                                            subindexStartforArea = indexMegaMatrixXaxis + 1;
                                        end
                                        FlagSetforHeader = 1;
                                        peakchange = peakchange + 1;
                                        % We want to write just one time Mass Hunter Peak Info e.g. (Peak RT Area Height Type Width FWHM) so its a flag
                                        if peakchange == 1
                                            indexMegaMatrixXaxis = indexMegaMatrixXaxis + 1;
                                            MegaMatrix(indexMegaMatrixXaxis, 3:5) = MassHunterDataFile(indexMassHunterFileData,1:3);
                                        end
                                        MegaMatrix(indexMegaMatrixXaxis, indexMegaMatrixYaxis:indexMegaMatrixYaxis+1) = [string(peak), ""];
                                        indexMegaMatrixYaxis = indexMegaMatrixYaxis + 2;
                                        subindexforArea = subindexforArea + 1;
                                    else % This is for Oxidized Peak Information  (MAKE SEPRATE FUNCTION FOR FORMATING)
                                        mod_pos=string(mascotfile(j,pep_var_mod_pos_column));
                                        size_mod_pos=length(regexpi(mod_pos{1}, '[0-9]'));

                                        ResidueInfo = [];
                                        %Extracting Oxidation information
                                        for indexOxiInfo=3:size_mod_pos
                                            if string(extractBetween(mod_pos,indexOxiInfo,indexOxiInfo))~='0'   %is entering the loop at 0. However it should enter when the char is 1
                                                pos=indexOxiInfo-2;
                                                ResidueInfo = [ResidueInfo,pos];

                                            end
                                        end
                                        if size(ResidueInfo,2) == 1
                                            %Finding Amino Acid with its position
                                            amino_acid = string(Sequence(ResidueInfo));
                                            %amino_acid=extractBetween(Sequence,pos,pos);
                                            position=string(ResidueInfo(1,1));
                                            aminoacidwithPos = strcat(amino_acid,position)
                                        else
                                            aminoacidwithPos = "";
                                            for iter=1: size(ResidueInfo,2)
                                                amino_acid = string(Sequence(ResidueInfo(1,iter)));
                                                position=string(ResidueInfo(1,iter));
                                                %We Want this M5+M12 But not his M5+M12+
                                                if iter ~= size(ResidueInfo,2)
                                                    aminoacidwithPos =strcat(aminoacidwithPos, amino_acid,position, '+');
                                                else
                                                    aminoacidwithPos =strcat(aminoacidwithPos, amino_acid,position);
                                                end
                                            end
                                        end

                                        if FlagSetforHeader == 0
                                            MegaMatrix(indexMegaMatrixXaxis, 1:2) = header;
                                        end
                                        FlagSetforHeader = 1;
                                        peakchange = peakchange + 1;
                                        % We want to write just one time Mass Hunter Peak Info e.g. (Peak RT Area Height Type Width FWHM) so its a flag
                                        if peakchange == 1

                                            indexMegaMatrixXaxis = indexMegaMatrixXaxis + 1;
                                            MegaMatrix(indexMegaMatrixXaxis, 3:5) = MassHunterDataFile(indexMassHunterFileData,1:3);
                                        end
                                        MegaMatrix(indexMegaMatrixXaxis, indexMegaMatrixYaxis:indexMegaMatrixYaxis+1) = [string(peak), aminoacidwithPos];
                                        indexMegaMatrixYaxis = indexMegaMatrixYaxis + 2;

                                    end
                                end
                            end
                        end
                    end
                    indexMegaMatrixXaxis = indexMegaMatrixXaxis + 1
                end
                ResultsAfterMatch=MegaMatrix;
                % Step 4: Calculate the Monoisotopic Theoretical mass of
                % Respective peptide
                %function to calculate the theoretical m/z , oxidize m/z and unoxidize m/z values
                [Oxidize_mz, Unoxidize_mz, Theoretical_peptide_weight] = MolecularWeight( Sequence )
                Oxidize_mz;
                Unoxidize_mz;
                % function to estimate the blocks of table that contains unoxidized
                % and oxidized mz information ( Peak, area, Retemtion time , matches)
                % Step 5 : sorting thta table Such that the Unoxidized Peak,
                % RT and Area comes 1nthe first  few rows Followed by
                % Oxidized Mz Data
                [Data_file_sort,File_AminoAcidHeader] =  Blocks(ResultsAfterMatch,Oxidize_mz,Unoxidize_mz)
                Data_file_sort;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % aligning the sequence of the peptide with the fasta file so that we can
                % get the start position of peptide sequence
                num=0
                pos=strfind(wholeSeq,string(Sequence))
                num=pos;
                rowNum= length(Data_file_sort)
                colNum= size(Data_file_sort)
                colnum=colNum(2)
                % where i is the iterative index to move along the sequence string
                for i= 1: length(Sequence)
                    data{i}=strcat(Sequence(i),num2str(i))

                    data_correct{i}=strcat(Sequence(i),num2str(i+num))

                    data_AminoAcid{2,colnum+i}=data{i}
                    data_AminoAcid{1,colnum+i}=data_correct{i}
                end
                %%%% data_aminoacid contains the header of the residues in that peptide
                data_AminoAcid
                %%%%%%%%%%%%%%%% merging 2 tables
                Data_file_sort

                if ~isfolder(strcat(pwd,'\\Results_matched'))
                    mkdir('Results_matched');
                end
                cd ('Results_matched')
                % Making directories as the name of peptide and sub
                % directories with the name of Replicate
                mkdir(CurrentPeptide)
                cd(CurrentPeptide)
                mkdir(CurrentReplicate)
                cd(CurrentReplicate)

                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],data_AminoAcid, 'Sheet1', 'C1');
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],Data_file_sort,'Sheet1','A3');
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],Data_file_sort,'Sheet1','B3');
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],Data_file_sort(:,5),'Sheet1','A3');

                FileProcessed=dir(fullfile(pwd,'\*.xlsx'))
                FileProcessedName = FileProcessed.name;
                % Extracting the complete name (Location + Name)
                FullFileName = fullfile(FileProcessed.folder, FileProcessedName);
                % Read the .xlsx file
                [num, txt, FileProcessed] = xlsread(FullFileName);

                cd ..\..\..
                % Step 6 : Populating the Matched Oxidized Residue table and Peak Split
                % Column
                FileProcessed=string(FileProcessed)
                % mock = string(mock);
                SizeOfFile=size(FileProcessed)
                SizeOfFile=SizeOfFile(2)
                [Index_var, colname] =  Columns(FileProcessed)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %penalty=1
                %%%%%%% assigning 1 under the oxidized residues
                %for  i= Index_var(1): Index_var(length(Index_var))
                % for  j= 3: size(FileProcessed)
                % for SeperateNumericals= 1:SizeOfFile
                % if FileProcessed(j,SeperateNumericals)== string(colname(i))
                %FileProcessed(j,i)= penalty
                %penalty=penalty+1
                % end
                % end
                %penalty=1
                %end
                % end
                %score=1
                %penalty=1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %fill mising value by zero
                cons="0"
                FileProcessed(3: size(FileProcessed),Index_var(1): Index_var(length(Index_var)))= fillmissing(FileProcessed(3: size(FileProcessed),Index_var(1): Index_var(length(Index_var))), 'constant',cons)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for  i= Index_var(1): Index_var(length(Index_var))
                    for  j= 3: size(FileProcessed)
                        for SeperateNumericals= 1:SizeOfFile
                            SplitTwoResidue=split(FileProcessed(j,SeperateNumericals),["+"])
                            for SplitResidueInd= 1:length(SplitTwoResidue)
                                if SplitTwoResidue(SplitResidueInd)==string(colname(i))

                                    FileProcessed(j,i)= str2double(FileProcessed(j,i)) + 1

                                end

                            end

                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% not sure
                alphabet=["A","B","C","D"]
                Alphaindex=1
                Headerr=FileProcessed(2,:)

                for  j= 3: size(FileProcessed)
                    Alphaindex=1
                    for SeperateNumericals= 1:SizeOfFile
                        if  contains(FileProcessed(j,SeperateNumericals),["+"])
                            SplitTwoResidue=split(FileProcessed(j,SeperateNumericals),["+"])
                            for SplitResidueInd= 1:length(SplitTwoResidue)
                                INDEX= find(contains(Headerr,SplitTwoResidue(SplitResidueInd)))

                                FileProcessed(j,INDEX)= FileProcessed(j,INDEX)+ alphabet(Alphaindex)

                            end
                            Alphaindex=Alphaindex+1

                        end

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%% assignining oxidation number in column 10 where \RowIndex is the row index
                totalrows= size(FileProcessed)
                totalrows= totalrows(1)
                for RowIndex= 1:size(FileProcessed)
                    if   contains(FileProcessed(RowIndex,6),string(Oxidize_mz(1:9)))
                        RowNumber= RowIndex+1
                        while ismissing(FileProcessed(RowNumber,6))
                            FileProcessed(RowNumber,10)=1;
                            RowNumber=RowNumber+1
                            if RowNumber >totalrows
                                break
                            end

                        end
                    else if contains(FileProcessed(RowIndex,6),string(Oxidize_mz(10:18)))
                            RowNumber= RowIndex+1
                            while ismissing(FileProcessed(RowNumber,6))
                                FileProcessed(RowNumber,10)=2;
                                RowNumber=RowNumber+1
                                if RowNumber > totalrows
                                    break
                                end

                            end
                    else if contains(FileProcessed(RowIndex,6),string(Oxidize_mz(19:27)))
                            RowNumber= RowIndex+1
                            while ismissing(FileProcessed(RowNumber,6))
                                FileProcessed(RowNumber,10)=3;
                                RowNumber=RowNumber+1
                                if RowNumber  > totalrows
                                    break
                                end

                            end
                    end
                    end
                    end
                end

                % connt the number of residues Oxidized in each peak
                for  i= 3:totalrows
                    new_dat=FileProcessed(i,Index_var(1):Index_var(length(Index_var)))
                    CountNumberOfResidues=sum(new_dat'~="0")
                    % IF the the total number of residues in a peak is
                    % greater than the oxidation number than place that
                    % number in a sepeprate column; this column will be
                    % identified as peak split
                    if  CountNumberOfResidues > str2double(FileProcessed(i,10))
                        FileProcessed(i,11)=FileProcessed(i,10);
                    end
                end

                %% Remove 'TC' or any other variable from 1st column and place '.x' ( TC0 --> 0.x)
                for i=3: totalrows
                    ColIdxTable=split(FileProcessed(i,1),[" "])
                    if ~ismissing(FileProcessed(i,1))
                        %If Row has a non missing value save that value in
                        % str variable
                        str = ColIdxTable(1);
                        SeperateNumericals = regexp(str,'(\d*)([a-z]*)','tokens');
                        OutputFileReplacedBy_X = [SeperateNumericals{:}];
                        OutputFileReplacedBy_X = OutputFileReplacedBy_X(~cellfun(@isempty,OutputFileReplacedBy_X));
                        OutputFileReplacedBy_X= append(OutputFileReplacedBy_X,{'.x'})
                        FileProcessed(i,1)=OutputFileReplacedBy_X
                    end
                end
                for i=3: totalrows
                    ColIdxTable=split(FileProcessed(i,2),[" "])
                    if ~ismissing(FileProcessed(i,2))
                        %If Row has a non missing value save that value in
                        % str variable
                        str = ColIdxTable(1);
                        SeperateNumericals = regexp(str,'(\d*)([a-z]*)','tokens');
                        OutputFileReplacedBy_X = [SeperateNumericals{:}];
                        OutputFileReplacedBy_X = OutputFileReplacedBy_X(~cellfun(@isempty,OutputFileReplacedBy_X));
                        OutputFileReplacedBy_X= append(OutputFileReplacedBy_X,{'.x'})
                        FileProcessed(i,1)=OutputFileReplacedBy_X
                    end
                end

                % remove the duplicate and empty rows

                [~,uidx] = unique(FileProcessed(:,1),'stable');
                FileWithoutDuplication = FileProcessed(uidx,:);
                DupliactionRemovedTable=FileWithoutDuplication
                totalrowsDUP= size(DupliactionRemovedTable)
                totalrowsDUP= totalrowsDUP(1)
                out=DupliactionRemovedTable(:,any(~ismissing(DupliactionRemovedTable(3:totalrowsDUP,:))))
                [Index_var, colname] =  Columns(DupliactionRemovedTable)

                % removing unnecessary Columns
                FinalData=DupliactionRemovedTable(:,[1,5,9,11,Index_var])
                FinalData(2,:) = []
                FinalData2= FinalData(any(~ismissing(FinalData),2),:)
                cd ('Results_matched')
                cd(CurrentPeptide)
                cd(CurrentReplicate)
                delete(FullFileName)
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],FinalData2, 'Sheet1', 'A1');
                cd ..\..\..
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This part reads the Intermediate file in .xlsx format and analyze it
% to get F value i.e. the proportion of the oxidized Residue and Area under
% the curve

% Step 1: Get .xlsx (Excel) files from the user. Note that these are the
% files that are the result of the module 1
% Specify the folder where the files live.
%myFolder='D:\DR. SHAHID\11March_complexity\Result\New folder\ResultsFromMatching'
myFolder= char(strcat(pwd,strcat('\','Results_matched')))
% Get a list of all files in the folder with the desired file name pattern.
FilePattern = fullfile(myFolder, '**\*.xlsx');
TheFiles = dir(FilePattern);
% Step 2: Extracting the File Name, Which is later used to make a Result folder of
% that particular peptide file
% for each peptide file present in the folder
if ~isfolder(strcat(pwd,'\\Result'))
    mkdir('Result');
end
for peptide = 1 : length(TheFiles)
    % Extracting the name of file with extension
    BaseFileName = TheFiles(peptide).name;
    % Extracting the complete name (Location + Name)
    FullFileName = fullfile(TheFiles(peptide).folder, BaseFileName);
    % Read the .xlsx file
    [num, txt, File] = xlsread(FullFileName);
    % Extracting the Name of file without the extension
    [myFolder,name,ext] = fileparts(FullFileName);
    fprintf(1, 'Now reading %s\n', FullFileName);
    % File=File(:,any(~ismissing(File(2:length(File),:))))
    %Get size of file
    SizeOfFile = size(File);
    % replicate names later used in naming the files
    rep= ["r1", "r2","r3"]
    %Convert file to a string to get headers in the file
    File = string(File);
    clear FileRemoveZeros
    clear ResidueOxidationfile
    %replace(File  (:,Index_var(1):Index_var(length(Index_var))) , "0",'')
    %F=replace(File  (:,:) , "0",'')
    [row_ind, row_end] = XRayDosageDataBlockIndices(File);

    % Step 4: Extract the  names of oxidized residues
    [Index_var, colname] =  NumberOfColumns(File)
    %col contains information about peak spliting
    col= File(:,4);
    col=str2double(col)
    col(isnan(col))=0;
    FileRemoveZeros(:,:) = regexprep(File(:,:), "0" ,'');
    empties = cellfun('isempty',FileRemoveZeros);
    FileRemoveZeros(empties) = {NaN}
    TotalRows=size(FileRemoveZeros)
    TotalRows=TotalRows(1)
    File=File(:,any(~ismissing(FileRemoveZeros(2:TotalRows,:))))
    % Step 3: Extract indices of data blocks at various X Ray dosages
    [Index_var, colname] =  NumberOfColumns(File)

    % Converting string file to doubble
    File_Main = str2double(File);

    %COVERTING NANS TO ZERO
    File_Main(isnan(File_Main))=0;

    %% Aminoacid and there reactivity value
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];
    reactivity= [ 0.14,2.9, 0.44 ,0.42,29.2,0.66,0.69,0.04,9.3,4.4, 4.4,2.2 ,20.5,11.2,1.0,1.4,1.6, 17.4,12.0,1.9]
    ThreeAmino_acids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"];
    % Matrix generated to save the caculation
    %Result1= zeros(6,17);  %for replicate 1
    % Mean and standard error of replicates error
    %Tab_Error= zeros(5,(length(Index_var)* 2))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each block
    for pos=1:length(row_ind)
        row_ind(pos);
        %Extracting rows between the index of first row and the last row within a block of X-ray dosage
        row= row_ind(pos): row_end(pos);
        length(row);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 6: Calculating of Sum of Areas of Each Replicate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sum of Areas of Replicate 1
        % file_main(row,5) contains spcific number of rows in column 2 - Unoxidized area(R1)
        format long
        total_1 = sum(File_Main(row,2));
        % file_main(row,15) contains  number of rows in column 3 - oxidized area(R1)
        format long
        total_2 = sum(File_Main(row,3));
        % sum of columns - oxidized area(R1) and Unoxidized area(R1)
        R1_total=total_1+ total_2;
        % sum is saved in 1st column of 'Result1' matrix
        Result1(pos,1)= R1_total;
        % Step 7: Calculating F values for each oxidized residue in a peptide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% for each Replicate   %%%%%%%%%%%%%%%%%%%%%%
        ResidueOxidationfile(:,Index_var)=File(:,Index_var)
        % Defining area as column 3 - oxidized area(R1)
        Area= File_Main(:,3);
        % Transpose of variable 'area' is taken for computataion purposes
        Area=Area';
        %declaring a variable count to calculate the total area of a spicified residue
        count=0;
        % reactivity value of 1st aminoacid
        RC_1a=0;
        % reactivity value of 2nd aminoacid
        RC_1b=0;
        RC_1c=0;
        name_col=[];
        Amino_acid1=0;
        Amino_acid2=0;
        Index=2;        % index variable  here is define as the number of column to store data

        % For each residue in the column
        for j=Index_var(1):Index_var(length(Index_var))

            % For each row in the selected area
            for i=row_ind(pos): row_end(pos)
                % checking wheter the residue has overlap retention time with
                % any other residues
                % checking wheter the residue has overlap retention time with
                % any other residues
                if col(i)==1 && ResidueOxidationfile(i,j)~= "0"
                    %Extracting the column name to Get the name of the amino acid
                    %and its reacticity.
                    name_col= convertStringsToChars(colname(j));
                    Amino_acid1=name_col(1);
                    IndexofResidue = strfind(amino_acids,Amino_acid1);
                    % Calculating the Ocuurance of that aminoacid in Peak
                    NumberOfOccurance=convertStringsToChars(ResidueOxidationfile(i,j))
                    NumberOfOccurance=str2double(NumberOfOccurance(1))
                    RC_1a= RC_1a + reactivity(IndexofResidue);
                    RC_1a=  NumberOfOccurance* RC_1a
                    % checking the neighbouring columns to identify the overlapped residue
                    %Fetching all the non zero numbers in the row under consideration
                    SingleRow= ResidueOxidationfile(i,:)~= "0" &  ~ismissing(ResidueOxidationfile(i,:))
                    TotalNonMissingValues=find(SingleRow~= 0 )
                    % If the  other non zeroalue is the same as index J tha do nothing
                    for TotalNonMissingValuesID= 1: length(TotalNonMissingValues)
                        if  TotalNonMissingValues(TotalNonMissingValuesID)== j
                            %do nothing
                            % if the index of other no zero value is other tha j
                            % than find the residue name its number of occurance
                        else
                            AminoID= TotalNonMissingValues(TotalNonMissingValuesID)
                            name_col= convertStringsToChars(colname(AminoID));
                            Amino_acid2=name_col(1);
                            IndexofResidue = strfind(amino_acids,Amino_acid2);

                            NumberOfOccurance=convertStringsToChars(ResidueOxidationfile(i, AminoID))
                            NumberOfOccurance=str2double( NumberOfOccurance(1))
                            RC_1b= RC_1b + reactivity(IndexofResidue);
                            RC_1b=  NumberOfOccurance* RC_1b
                        end
                    end
                    % sum of reactivitity constants of 2 aminoacids
                    total_RC= RC_1a+ RC_1b;
                    RC=RC_1a/total_RC;
                    fprintf('Reactivity of %s with %s at row %f is %f .\n ',Amino_acid1,Amino_acid2, i,RC)
                    % Divided area of residues having same Retention time
                    Total_area_R=Area(i)* RC;
                    % Add the area to the variable count
                    count=count+Total_area_R;
                    % Sum of area of  residue that doesnot operlap with other residues
                    % checking wheter the residue has overlap retention time with
                    % any other residues
                    % checking wheter the residue has overlap retention time with
                    % any other residues
                    % Extarctng the Area If the peak has oxidation number 2
                else if col(i)==2 && ResidueOxidationfile(i,j)~= "0"
                        %Extracting the column name to Get the name of the amino acid
                        %and its reacticity.
                        name_col= convertStringsToChars(colname(j));
                        Amino_acid1=name_col(1);
                        IndexofResidue = strfind(amino_acids,Amino_acid1);
                        RC_1a= RC_1a + reactivity(IndexofResidue);
                        NumberOfOccurance_a=convertStringsToChars(ResidueOxidationfile(i,j))
                        NumberOfOccurance_a=str2double(NumberOfOccurance_a(1))
                        %estimating the different residues of overlapped peak
                        % calculating how many times this residue is present
                        PositionOfResidue=convertStringsToChars(ResidueOxidationfile(i,j))
                        length(PositionOfResidue)
                        for len= 2: length(PositionOfResidue)
                            UniqueAminoAcidString(len-1)= PositionOfResidue(len)

                        end
                        % checking the neighbouring columns to identify the overlapped residue
                        %Fetching all the non zero numbers in the row under consideration
                        SingleRow= ResidueOxidationfile(i,:)~= "0" &  ~ismissing(ResidueOxidationfile(i,:))
                        TotalNonMissingValues=find(SingleRow~= 0 )
                        % checking if the nonzero number index is same as index
                        % 'j'
                        for TotalNonMissingValuesID= 1: length(TotalNonMissingValues)
                            if  TotalNonMissingValues(TotalNonMissingValuesID)== j
                                % do nothing
                            else
                                %if the nonzero index is not same as "j" than
                                %find the names and accurances of other
                                %residues in that peak
                                AminoID= TotalNonMissingValues(TotalNonMissingValuesID)
                                name_col= convertStringsToChars(colname(AminoID));
                                Amino_acid2=name_col(1);
                                AminoAcideResidueTwo(TotalNonMissingValuesID)=Amino_acid2
                                IndexofResidue = strfind(amino_acids,Amino_acid2);
                                RC_1b(TotalNonMissingValuesID)= reactivity(IndexofResidue);
                                %% getting the occurances of the residues
                                PositionOfResidue=convertStringsToChars(ResidueOxidationfile(i,AminoID))
                                length(PositionOfResidue)
                                for len= 2: length(PositionOfResidue)
                                    UniqueAminoAcidString=append(UniqueAminoAcidString,PositionOfResidue(len))
                                end
                            end
                        end
                        % A unique aray that contains all the possible numbers
                        % a residue can be present
                        UniqueAminoAcidString=unique(UniqueAminoAcidString)
                        Idrow=1

                        for lengthofCheck=1: length(UniqueAminoAcidString)
                            if contains(ResidueOxidationfile(i,j),UniqueAminoAcidString(lengthofCheck))
                                denom(Idrow,lengthofCheck)= RC_1a
                            end
                        end
                        % this block finds the Reactivities of the amino acids and place them in
                        % matrix 'denom' such that each column contains the residues whoes
                        % reactivities are to multiplied
                        for TotalNonMissingValuesID= 1: length(TotalNonMissingValues)
                            if  TotalNonMissingValues(TotalNonMissingValuesID)== j
                            else
                                Idrow=Idrow+1
                                AminoID= TotalNonMissingValues(TotalNonMissingValuesID)
                                name_col= convertStringsToChars(colname(AminoID));
                                Amino_acid2=name_col(1);
                                AminoAcideResidueTwo(TotalNonMissingValuesID)=Amino_acid2
                                IndexofResidue = strfind(amino_acids,Amino_acid2);
                                RC_1b= reactivity(IndexofResidue);
                                for lengthofCheck=1: length(UniqueAminoAcidString)
                                    if contains(ResidueOxidationfile(i,AminoID),UniqueAminoAcidString(lengthofCheck))
                                        denom(Idrow,lengthofCheck)= RC_1b
                                    end
                                end
                            end
                        end
                        TableForDenominator = denom;
                        % For the even-numbered rows, circularly shift the elements
                        TableForDenominator(~TableForDenominator)=1  % replace 0 elements by 1
                        Product_output=prod(TableForDenominator,1)
                        Denominator=0
                        % calculating neumerator
                        Numerator= RC_1a*NumberOfOccurance_a

                        % calculating denominator
                        for productId=1:length(Product_output)
                            Denominator=  Denominator+ Product_output(productId)
                        end
                        RC=Numerator/Denominator;
                        % Divided area of residues having same Retention time
                        Total_area_R=Area(i)* RC;
                        % Add the area to the variable count
                        count=count+Total_area_R;
                        %%%%% if peak split has variable 3
                else if col(i)==3 && ResidueOxidationfile(i,j)~= "0"
                        %Extracting the column name to Get the name of the amino acid
                        %and its reacticity.
                        name_col= convertStringsToChars(colname(j));
                        Amino_acid1=name_col(1);
                        IndexofResidue = strfind(amino_acids,Amino_acid1);
                        RC_1a= RC_1a + reactivity(IndexofResidue);
                        NumberOfOccurance_a=convertStringsToChars(ResidueOxidationfile(i,j))
                        NumberOfOccurance_a=str2double(NumberOfOccurance_a(1))
                        %estimating the different residues of overlapped peak
                        % calculating how many times this residue is present
                        PositionOfResidue=convertStringsToChars(ResidueOxidationfile(i,j))
                        length(PositionOfResidue)

                        for len= 2: length(PositionOfResidue)
                            UniqueAminoAcidString(len-1)= PositionOfResidue(len)
                        end
                        % checking the neighbouring columns to identify the overlapped residue
                        %Fetching all the non zero numbers in the row under consideration
                        SingleRow= ResidueOxidationfile(i,:)~= "0" &  ~ismissing(ResidueOxidationfile(i,:))
                        TotalNonMissingValues=find(SingleRow~= 0 )
                        % checking if the nonzero number index is same as index
                        % 'j'
                        for TotalNonMissingValuesID= 1: length(TotalNonMissingValues)
                            if  TotalNonMissingValues(TotalNonMissingValuesID)== j
                                % do nothing
                            else
                                %if the nonzero index is not same as "j" than
                                %find the names and accurances of other
                                %residues in that peak
                                AminoID= TotalNonMissingValues(TotalNonMissingValuesID)
                                name_col= convertStringsToChars(colname(AminoID));
                                Amino_acid2=name_col(1);
                                AminoAcideResidueTwo(TotalNonMissingValuesID)=Amino_acid2
                                IndexofResidue = strfind(amino_acids,Amino_acid2);
                                RC_1b(TotalNonMissingValuesID)= reactivity(IndexofResidue);
                                %% getting the occurances of the residues
                                PositionOfResidue=convertStringsToChars(ResidueOxidationfile(i,AminoID))
                                length(PositionOfResidue)
                                for len= 2: length(PositionOfResidue)
                                    UniqueAminoAcidString=append(UniqueAminoAcidString,PositionOfResidue(len))
                                end
                            end
                        end
                        % A unique aray that contains all the possible numbers
                        % a residue can be present
                        UniqueAminoAcidString=unique(UniqueAminoAcidString)
                        Idrow=1
                        for lengthofCheck=1: length(UniqueAminoAcidString)
                            if contains(ResidueOxidationfile(i,j),UniqueAminoAcidString(lengthofCheck))
                                denom(Idrow,lengthofCheck)= RC_1a
                            end
                        end
                        % this block finds the Reactivities of the amino acids and place them in
                        % matrix 'denom' such that each column contains the residues whoes
                        % reactivities are to multiplied

                        for TotalNonMissingValuesID= 1: length(TotalNonMissingValues)
                            if  TotalNonMissingValues(TotalNonMissingValuesID)== j
                            else
                                Idrow=Idrow+1
                                AminoID= TotalNonMissingValues(TotalNonMissingValuesID)
                                name_col= convertStringsToChars(colname(AminoID));
                                Amino_acid2=name_col(1);
                                AminoAcideResidueTwo(TotalNonMissingValuesID)=Amino_acid2
                                IndexofResidue = strfind(amino_acids,Amino_acid2);
                                RC_1b= reactivity(IndexofResidue);
                                for lengthofCheck=1: length(UniqueAminoAcidString)
                                    if contains(ResidueOxidationfile(i,AminoID),UniqueAminoAcidString(lengthofCheck))
                                        denom(Idrow,lengthofCheck)= RC_1b
                                    end
                                end
                            end
                        end
                        TableForDenominator = denom;
                        % For the even-numbered rows, circularly shift the elements
                        TableForDenominator(~TableForDenominator)=1  % replace 0 elements by 1
                        Product_output=prod(TableForDenominator,1)
                        Denominator=0
                        % calculating neumerator
                        Numerator= RC_1a*NumberOfOccurance_a
                        % calculating denominator
                        for productId=1:length(Product_output)
                            Denominator=  Denominator+ Product_output(productId)
                        end
                        RC=Numerator/Denominator;
                        % Divided area of residues having same Retention time
                        Total_area_R=Area(i)* RC;
                        % Add the area to the variable count
                        count=count+Total_area_R;
                        %If peak split doesnot have any information
                else if  col(i)==0 && ResidueOxidationfile(i,j)~= "0"
                        SingleRow= ResidueOxidationfile(i,:)~= "0" &  ~ismissing(ResidueOxidationfile(i,:))
                        TotalNonMissingValues=find(SingleRow~= 0 )

                        NumberOfOccurance_a=convertStringsToChars(ResidueOxidationfile(i,j))
                        NumberOfOccurance_a=str2double(NumberOfOccurance_a(1))
                        count=count+(Area(i)* NumberOfOccurance_a);

                end
                end
                end
                end

                RC_1a=0;
                RC_1b=0;
                RC_1c=0;
                clear  denom

            end
            fprintf('Value of R curve of %s is %f .\n ', Amino_acid1,count)
            % F_1 is the ratio of oxidizes residues
            f_1= count/Result1(pos,1);
            % final_1 is the reverse of F_1 that is (1_F)
            final_1= 1 - f_1 ;
            % Storing the result of counts for each residues in matrix Result
            Result1(pos,Index)= count;
            Index=Index+1;
            %Storing the result of oxidize Fraction for each residues in matrix Result
            Result1(pos,Index)= f_1;
            Index=Index+1;
            % Storing the result of 1-F for each residues in matrix Result
            Result1(pos,Index)= final_1;
            Index=Index+1;
            count=0;
            RC_1a=0;
            RC_1b=0;

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 8: generating excel files of result matrices with proper Header
    % Result1 matrix contain calculations for replicate 1

    col_header=[];
    idxCol = 0;
    idxNew = 1;
    % Defining Column Header by using the variable 'colname' extracted in
    % step 4. This variable contain names of the oxidized residues
    for idx = 1:size(colname,2)
        idxCol = idxCol + 1;
        if ~ismissing(colname(1,idxCol))
            col_header{idxNew} =  colname (1,idxCol) + '_SUM';
            col_header{idxNew+1} = colname (1,idxCol) + '_F';
            col_header{idxNew+2} = colname (1,idxCol) + '_1-F';
            idxNew = idxNew + 3;
        end
    end
    % Defining Row Header
    row_header= { 'Dose'; '0.x';'10.x';'25.x';'50.x';'75.x';'100.x'};
    %FileName='Result_'+ rep(peptide);
    %save(Filename,result1)
    cd ('Result')
    %  Creating a folder with the name of the peptide file name extracting from
    % the Step 2
    BaseFileName1=BaseFileName(1:8)

    if ~isfolder(strcat(pwd, BaseFileName1))
        mkdir(BaseFileName1);
    end
    cd (BaseFileName1)
    if ~isfolder(strcat(pwd,name))
        mkdir(name);
    end
    % Change the working directory to the newly made directory inorder to store
    % results
    cd (name)
    % Writing the matrix created in Step 7 in excel file.
    xlswrite('Replicate.xlsx',Result1,'Sheet1','B2');
    xlswrite('Replicate.xlsx',col_header,'Sheet1','C1');
    xlswrite('Replicate.xlsx',row_header,'Sheet1','A1');

    cd ..\..\..
    clear Result1
end
% Reading the Newly made result Directory to Read the excel files in the
% non empty subfolders
cd ('Result')
ResultDirectory = dir(pwd);
NonEmptyFolderlength= sum([ResultDirectory(~ismember({ResultDirectory.name},{'.','..'})).isdir])
for IndexofFolder=3:NonEmptyFolderlength+2
    Subfolder=ResultDirectory(IndexofFolder)
    FullFileNameofPeptide = fullfile(ResultDirectory(IndexofFolder).folder, Subfolder.name)
    cd (FullFileNameofPeptide)
    rootdir = pwd

    %get list of files and folders in any subfolder
    TheFiles = dir(fullfile(rootdir, '**\*.xlsx'));
    TheFiles = TheFiles(~[TheFiles.isdir])
    for peptide = 1 : length(TheFiles)
        % Extracting the name of file with extension
        BaseFileName = TheFiles(peptide).name;
        FullFileName = fullfile(TheFiles(peptide).folder, BaseFileName);
        % Read the .xlsx file
        [num, txt,ArrayCellofReplicates{peptide}] = xlsread(fullfile(FullFileName))
        %Extracting the Name of file without the extension
        [myFolder,Name,ext] = fileparts( BaseFileName);
        fprintf(1, 'Now reading %s\n',  BaseFileName);
    end


    % Reading the Results of 3 rweplicates of a single peptide into tables
    % with different names
    ResultRep1=ArrayCellofReplicates{1,1};
    ResultRep1=string(ResultRep1)
    ResultRep1=str2double(ResultRep1)
    ResultRep2=ArrayCellofReplicates{1,2};
    ResultRep2=string(ResultRep2)
    ResultRep2=str2double(ResultRep2)
    ResultRep3=ArrayCellofReplicates{1,3};
    ResultRep3=string(ResultRep3)
    ResultRep3=str2double(ResultRep3)
    [rownum,colnum]=size(ResultRep1)
    [rownum1,colnum2]=size(ResultRep2)
    newCol=[0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1]
    if colnum== colnum2
        % do nothing
    else
        ResultRep3(2:7,colnum2+1:colnum2+3)=newCol

        ResultRep2(2:7,colnum2+1:colnum2+3)=newCol
    end

    % Step 9: calculating the mean and Std error of 1-F value for each Residue of each peptide

    for row_2= 2:length(row_ind)+1
        ind_2=1;
        % for each column having value of 1_F in Result matrix
        for r_col= 5:3:length(ResultRep2)
            % calculating Mean and Standard error of a '1-F 'value of a  residue in both result
            % matrix
            Mean = mean([ResultRep1(row_2,r_col) ResultRep2(row_2,r_col) ResultRep3(row_2,r_col)], 'all');
            std_error= std2([ResultRep1(row_2,r_col) ResultRep2(row_2,r_col) ResultRep3(row_2,r_col)])/sqrt(3)
            rev= 1/std_error
            % Saving the result in Tab_error Matrix
            Tab_Error(row_2,ind_2)= Mean;
            ind_2=ind_2+1;
            Tab_Error(row_2,ind_2)=std_error;
            ind_2=ind_2+1;
            Tab_Error(row_2,ind_2)=rev;
            ind_2=ind_2+1;
        end
    end
    % Defining Column Header by using the variable 'colname' extracted in
    % step 4. This variable contain names of the oxidized residues

    head=ArrayCellofReplicates{1,1}
    TotalColumn= size(head)
    TotalColumn=TotalColumn(2)
    head= head(1, 3:TotalColumn)
    r = strtok(head, '_')
    col_names=unique(r,'stable')
    col_names=string(col_names)
    cd ..\..
    cd (FullFileNameofPeptide)
    header = [];
    idxCol = 0;
    idxNew = 1;
    for idx = 1:size(col_names,2)
        idxCol = idxCol + 1;
        if ~ismissing(col_names(1,idxCol))
            header{idxNew} =  col_names (1,idxCol) + '_MEAN';
            header{idxNew+1} = col_names (1,idxCol) + '_ERROR',
            header{idxNew+2} = col_names (1,idxCol) + '_rev_error';
            idxNew = idxNew + 3;
        end
    end
    % generating excel files of result matrices with proper Header
    %tab_error Matrix contain calculations for Mean and Error
    xlswrite('Mean_R1_R2_R3',Tab_Error,'Sheet1','B1');
    xlswrite('Mean_R1_R2_R3',header,'Sheet1','B1');
    xlswrite('Mean_R1_R2_R3',row_header,'Sheet1','A1');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step read the result file in xls format that we have saved in  step 9
    FileMean=readtable('Mean_R1_R2_R3.xls');
    % Convert it to table for easy indexing

    FileMean = table2cell(FileMean);
    f = @(x) any(isnan(x));
    B = cellfun(f, FileMean);
    C = all(B);
    FileMean(:,C) = []
    Dosage= [0;10;25;50;75;100]
    cd ..\..\
    % Step 10: Fitting the exponential curve for the mean values
    [ExponentialFit, DataSlope] =  InitialSlopeMean( FullFileNameofPeptide,header,FileMean,File,col_names)
    % row names of the table comtaining slope values
    row_head= {'coefficients';'Lower';'Higher';'Average';'':'Gradient'};
    % Extracting column names from the function and writing the Slope values
    % in an excel file
    header = [];
    idxCol = 0;
    idxNew = 1;
    for idx = 1:size(col_names,2)
        idxCol = idxCol + 1;
        if ~ismissing(col_names(1,idxCol))
            header{idxNew} =  col_names (1,idxCol) + '_a';
            header{idxNew+1} = col_names (1,idxCol) + '_b',

            idxNew = idxNew + 2;
        end
    end
    DataSlope;
    cd ('Result')
    cd ( FullFileNameofPeptide)
    xlswrite('Slope',DataSlope,'Sheet1','B2');
    xlswrite('Slope',header,'Sheet1','B1');
    xlswrite('Slope',row_head,'Sheet1','A1');
    AA_header= col_names
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 11: identifying the oxidized Residues Intrinsic Reactivity and
    % Calculating Protection factor   ( PF= Intrinsic reactivity / Gradient)
    AminoAcidPosition=2;
    Row_number=1;

    for res=1:length(col_names)
        Name_AA=convertStringsToChars(AA_header(res));
        Amino_acidn=Name_AA(1);
        IndexofAminoacid = strfind(amino_acids,Amino_acidn);
        itr=reactivity(IndexofAminoacid);
        for Row_number= 5
            ProtectionFactor= itr / DataSlope(Row_number,AminoAcidPosition);
            Data_ProtectionFactor(2,AminoAcidPosition)=ProtectionFactor;
        end
        AminoAcidPosition=AminoAcidPosition+2;
    end
    row_head= {'protection_factor';'';'Average'};
    % Extracting column names from the function and writing the Slope values
    % in an excel file
    cd ../..
    % Generating Header
    header = [];
    idxCol = 0;
    idxNew = 1;
    for idx = 1:size(col_names,2)
        idxCol = idxCol + 1;
        if ~ismissing(col_names(1,idxCol))
            header{idxNew+1} = col_names (1,idxCol) + '_PF',

            idxNew = idxNew + 2;
        end
    end
    cd ('Result')
    cd ( FullFileNameofPeptide)
    xlswrite('Slope',Data_ProtectionFactor,'Sheet1','B10');
    xlswrite('Slope',header,'Sheet1','B9');
    xlswrite('Slope',row_head,'Sheet1','A9');
    cd ../..
    % Convert it to table for easy indexing
    SASA=cell(1, length(DataSlope));
    SASA(1,:) = {0};
    AA_colId=2;
    % to extract the AminoAcid resuide  and SASA value
    for res=1:length(col_names)
        Name_AA=convertStringsToChars(AA_header(res));
        Amino_acidname=Name_AA(1);
        Amino_acidNum=Name_AA(2:length(Name_AA));
        Amino_acidNum=str2num(Amino_acidNum);
        IndexofAminoacid = strfind(amino_acids,Amino_acidname);
        AA_symbol=ThreeAmino_acids(IndexofAminoacid);
        SASA(1,AA_colId) =FileSASA(Amino_acidNum-1,7);
        AA_colId=AA_colId+2;
    end
    title= {'SASA'};
    % Generating Header
    headerSASA = [];
    idxCol = 0;
    idxNew = 1;
    for idx = 1:size(col_names,2)
        idxCol = idxCol + 1;
        if ~ismissing(col_names(1,idxCol))
            headerSASA{idxNew+1} = col_names (1,idxCol) + ' ',

            idxNew = idxNew + 2;
        end
    end
    cd ('Result')
    cd ( FullFileNameofPeptide)
    xlswrite('Slope',SASA,'Sheet1','B20');
    xlswrite('Slope',title,'Sheet1','A19');
    xlswrite('Slope',headerSASA,'Sheet1','B19');
    %%%%%%%% Taking log10( -PF) of data
    logData=log10(-Data_ProtectionFactor)
    % defining header
    Column_IDX=1
    bb=5;
    for ColIdxTable=2:2:length(logData)
        PF_log(:,Column_IDX)=logData(:,ColIdxTable)
        SASA_X(1,Column_IDX)=SASA(:,ColIdxTable)
        Column_IDX=Column_IDX+1;
    end
    SASA_X=cell2mat(SASA_X)
    %% Sorting the data into X and Y axis for fit
    yaxis=PF_log(2,:);
    yaxis=yaxis',
    xaxis= SASA_X;
    xaxis=xaxis';
    DATA_FIT(:,1)=xaxis
    DATA_FIT(:,2)=yaxis
    DATA_FIT=sortrows(DATA_FIT)
    xaxis= DATA_FIT(:,1)
    yaxis= DATA_FIT(:,2)
    TotalCountRows=size(xaxis)
    TotalCountRows=TotalCountRows(1)
    logData=logData(2,:)

    title= {'LogPF'};
    xlswrite('Slope',logData,'Sheet1','B25');
    xlswrite('Slope',title,'Sheet1','A24');
    if TotalCountRows > 1

        % Linear Fit
        fit_Curve=fit(xaxis,yaxis,'poly1');


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Identifying 95% confidence In tervals for coefficients of Fit
        % int=confint(fit_Curve)
        %int(3,:)=coeffvalues(fit_Curve)

        figure (1)
        set(gcf,'Visible','off')
        % plotting the curve
        % plotting confidence interval in graph
        %y_fit = int(3,1)*exp(int(3,2)*xaxis);
        %y_cbl = int(1,1)*exp(int(1,2)*xaxis);
        %y_cbh = int(2,1)*exp(int(2,2)*xaxis);
        % plot(xaxis,y_fit,'r',xaxis,y_cbl,'b:',xaxis,y_cbh,'b:')
        %plot(ffff )
        plot(fit_Curve)
        hold on
        plot(xaxis ,yaxis,'c.','MarkerSize',15)
        ylim([0 6.5])
        hold off
        title_string = strcat('Log(PF) vs SASA');
        %title (title_string, 'Fontsize', 20 , 'fontweight', 'bold', 'Color', [0.8 0 0]);
        % title (title_string)
        xlabel('SASA', 'Fontsize', 16, 'fontweight', 'bold', 'Color', [0 0 0]);
        ylabel('Log(PF)', 'Fontsize', 16, 'fontweight', 'bold', 'Color', [0 0 0]);
        legend('Fit', 'Data')
        plotname = strcat('SASA.png');
        saveas(gcf,plotname)
    end
    cd ../../

    % To clear the variables
    clear SASA
    clear SASA_X
    clear data_slope
    clear Tab_Error
    clear logData
    clear Data_ProtectionFactor
    clear  SASA_X
    clear PF_log
    clear DATA_FIT



end
cd ('Result')
ResultDirectory = dir(pwd);
%number of non empty folders in directory or length of list in directory
NonEmptyFolderlength= sum([ResultDirectory(~ismember({ResultDirectory.name},{'.','..'})).isdir])
i=2;
for IndexofFolder=3:NonEmptyFolderlength+2
    Subfolder=ResultDirectory(IndexofFolder)
    % to Get the full name and location of Nonempty Directory
    FullFileNameofPeptide = fullfile(ResultDirectory(IndexofFolder).folder, Subfolder.name)
    cd (FullFileNameofPeptide)

    % Fetching and sorting the PF data and SASA data of Each peptide into a
    % single file which is alter used to compute the final output  graph and
    % Table
    rootdir = pwd
    %get list of files and folders in any subfolder
    TheFiles = dir(fullfile(rootdir, 'Slope.xls'));
    TheFiles = TheFiles(~[TheFiles.isdir])

    % Extracting the name of file with extension
    BaseFileName = TheFiles.name;
    FullFileName = fullfile(TheFiles.folder, BaseFileName);

    % Read the .xlsx file
    [num, txt,File] = xlsread(fullfile(FullFileName))
    row_ind = find(cellfun('length',regexp(string(File),'LogPF')) == 1);

    % convert file to string
    file = string(File)
    startIdx = row_ind(1)
    sizefile=size(File)
    SizeOfFile=sizefile(2)

    % find the missing values
    MissingTable = ismissing(file(startIdx:length(file)))
    idx = find(MissingTable~=0, 1, 'first')
    endID= (idx + startIdx) - 1
    % placing the output from SASA table and LogPF
    blocks(3,i: (i+SizeOfFile)-2)= File(endID,2:SizeOfFile)
    row_ind2 = find(cellfun('length',regexp(string(File),'SASA')) == 1);
    file = string(File)
    startIdx2 = row_ind2(1)
    sizefile=size(File)
    SizeOfFile=sizefile(2)
    MissingData = ismissing(file(startIdx2:length(file)))
    idx2 = find(MissingData~=0, 1, 'first')
    endID2= (idx2 + startIdx2) - 1
    blocks(1:2,i: (i+SizeOfFile)-2)= File( startIdx2:endID2,2:SizeOfFile)

    i= size(blocks,2)+1;
end
Tab=blocks(:,2:length(blocks))

out= string(Tab)
RemoveMissing=ismissing(out(1,:))
Tab=Tab(:,~RemoveMissing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sorting the input for Linear fit
yaxis=Tab(3,:);
yaxis=yaxis',
xaxis= Tab(2,:);
% Taking Transpose
xaxis=xaxis';
xaxis=cell2mat(xaxis)
yaxis=cell2mat(yaxis)
DATA_FIT(:,1)=xaxis
DATA_FIT(:,2)=yaxis
DATA_FIT=sortrows(DATA_FIT)
xaxis= DATA_FIT(:,1)
yaxis= DATA_FIT(:,2)

%Fitting linear regression to final conplied data of all the peptides and
%their Identified Oxidized Residues
Final_Linear_fit=fit(xaxis,yaxis,'poly1');
%CALCULATING 95% Prediction bounds
PredictionInterval =predint(Final_Linear_fit,xaxis,0.95,'functional','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identifying 95% confidence Intervals for coefficients of Fit
int=confint(Final_Linear_fit)
int(3,:)=coeffvalues(Final_Linear_fit)
figure (1)
set(gcf,'Visible','off')
% plotting the curve
cd ..
plot(Final_Linear_fit)
hold on
xlim([0 inf])
plot(xaxis ,yaxis,'c.','MarkerSize',15)
plot(xaxis,PredictionInterval,'m--')
hold off
title_string = strcat('Log(PF) vs SASA');
xlabel('SASA', 'Fontsize', 16, 'fontweight', 'bold', 'Color', [0 0 0]);
ylabel('Log(PF)', 'Fontsize', 16, 'fontweight', 'bold', 'Color', [0 0 0]);
legend({'Fit','Data','Lower Interval', 'Higher Interval'},'FontSize',4)
plotname = strcat('SASAmain.png');
saveas(gcf,plotname)
xlswrite('Slope_main',DATA_FIT,'Sheet1','B1');

% Split the Residue name into its alphabetic and numerical components( F8
% as  F and 8)
Gradient_SASA_Table=(Tab)'
for i=1: length(Gradient_SASA_Table)
    SplitResidueName=split(Gradient_SASA_Table(i,1),[" "])
    Gradient_SASA_Table(i,4)=regexprep(Gradient_SASA_Table(i,1),'\D','')
end
%%%%%%%% Removing un-necessary rows
Tableleft=cell2table(FileSASA(1:length(FileSASA),1:2), 'VariableNames', {'Position', 'Residue'})
Tableright=cell2table(Gradient_SASA_Table, 'VariableNames', {'Residue Symbol', 'SASA','LogPF','Position'})
TableOutput = convertvars(Tableleft,{'Position'},'string')
TableSASASort = convertvars(Tableright,{'Position'},'string')
Tableright.Position=str2double(Tableright.Position)
TableMerge=outerjoin(Tableleft,Tableright)
indx_pos=[1,2,4,5]
FinalTableOutput = TableMerge(:, indx_pos)
% write the output table
writetable(FinalTableOutput,'PF_SASA_tab.xls','Sheet',1);