%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM-XFMS: A MATLAB Toolbox to Analyze X-ray Footprinting Mass Spectrometry Data %
%                                    Version 1.0.0                                     %
%        Copyright (c) Biomedical Informatics & Engineering Research Laboratory,       %
%           Lahore University of Management Sciences Lahore (LUMS), Pakistan.          %
%                           (http://biolabs.lums.edu.pk/BIRL)                          %
%                                (safee.ullah@gmail.com)                               %
%                            Last Modified on: 26-July-2022                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FileMerge()
% This function reads the intermediate files and merge them on basis of the same header
myFolder= char(strcat(pwd,strcat('\','Results_matched')));
% Get a list of all files in the folder with the desired file name pattern.
FilePattern = fullfile(myFolder, '**\*.xlsx');
TheFiles = dir(FilePattern);
location=cell(10, 1);
for peptide = 1 : length(TheFiles)
    % Extracting the name of file with extension
    BaseFileName = TheFiles(peptide).name;
    % Extracting the complete name (Location + Name)
    FullFileName = fullfile(TheFiles(peptide).folder, BaseFileName);
    location(peptide)={fullfile(TheFiles(peptide).folder, BaseFileName)};
    % Read the .xlsx file;
    [num, txt, File{peptide}] = xlsread(FullFileName);
    % Extracting the Name of file without the extension
    [myFolder,name,ext] = fileparts(FullFileName);
    fprintf(1, 'Now reading %s\n', FullFileName);
    % File=File(:,any(~ismissing(File(2:length(File),:))))
    %Get size of file
    SizeOfFile = size(File);
end
length(File);
idxx=1;
%%%%% For replicate 1
for RowIndexMergeFile= 1:3: length(File)
    for ColumnIndexMergeFile= 1:3:length(File)
        [num, txt, FileMat] = xlsread(string(location(RowIndexMergeFile)));
        if RowIndexMergeFile== ColumnIndexMergeFile
            %do nothing

        elseif ColumnIndexMergeFile>RowIndexMergeFile
            [num, txt, FileMat2] = xlsread(string(location(ColumnIndexMergeFile)));
            %matching files
            logicalArray=ismember(string(FileMat2(1,:)),string(FileMat(1,:)));
            num=0;
            num2=0;
            % if Header is same
            if    any (logicalArray(:),1)
                % find  the index of blocks of 2 files
                [row_ind1, row_end1] = XRayDosageDataBlockIndicesCaseOne(FileMat);
                [row_ind2, row_end2] = XRayDosageDataBlockIndicesCaseTwo(FileMat2);

                for pos=1:length(row_ind1)
                    %Extracting rows between the index of first row and the last row within a block of X-ray dosage
                    row= row_ind1(pos): row_end1(pos);
                    length(row);
                    %%%%%%%%%%%%%%%%%%%%%%%
                    %saving the rows and columns having same dosage and
                    %same residue in new file
                    row2= row_ind2(pos): row_end2(pos);
                    index=5;
                    data(row_ind1(pos)+num: row_end1(pos)+num,1:4)= FileMat(row_ind1(pos): row_end1(pos),1:4);
                    data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,1:4)= FileMat2(row_ind2(pos): row_end2(pos),1:4);
                    FileMat2(row_ind2(pos): row_end2(pos),1:4)= {0};
                    % add the length
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% for residues
                        % it replaces Zeros at postiosn thta have beem merged
                    % so the ultimately the whole file has nothing but zero
                    %  and it cannot be resued again
                    for index2= 5: size(FileMat2,2)
                        if string(FileMat2(1,index2))== string(FileMat(1,index))
                            data(row_ind1(pos)+num: row_end1(pos)+num,index2)= FileMat(row_ind1(pos): row_end1(pos),index);
                            data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,index2)= FileMat2(row_ind2(pos): row_end2(pos),index2);
                            FileMat2(row_ind2(pos): row_end2(pos),index2)={0};
                            index=min(size(FileMat,2),index+1);
                        else
                            data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,index2)= FileMat2(row_ind2(pos): row_end2(pos),index2);
                            FileMat2(row_ind2(pos): row_end2(pos),index2)={0};
                        end
                    end
                    num=num+length(row2);
                    num2=num2+length(row);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %emptythe file 
                          FirstColumn = data(:,1);
                Storage = [];
                for index = 1:size(FirstColumn,1)
                    CellTemp = FirstColumn(index);
                    if iscellstr(CellTemp)
                        strTemp = string(CellTemp);
                        Storage = [Storage; [strTemp index] ];
                    end
                end
                idx=[];
                jidx=1;
                for id= 2:2:size(Storage,1)
                    idx(jidx)= Storage(id,2);
                    jidx=jidx+1;
                end
                data(idx,:)=[];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Make new folder to save newly made merged file
                if ~isfolder(strcat(pwd,'\\Resultsnew1'))
                    mkdir('Resultsnew1');
                end

                cd ('Resultsnew1')
                Header = string(FileMat2(1,:));
                %Location of the files that are merged 
                loc= location(RowIndexMergeFile);
                out=regexp(loc,'\','split');
                out=out{1};
                length(out);
                % creating the folder with appropriate peptide name
                CurrentPeptide=out{length(out)-2};
                CurrentReplicate=out{length(out)-1};
                mkdir(CurrentPeptide);
                cd(CurrentPeptide);
                mkdir(CurrentReplicate);
                cd(CurrentReplicate);
                directories{idxx}=CurrentPeptide;
                idxx=idxx+1;
                delete([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx']);
                temp = data(:,5: size(data,2));
                temp(cellfun(@(x)isscalar(x) && isnan(x),temp)) = {0};
                data(:,5: size(data,2)) = temp;
                tf=cellfun('isempty',data );
                data(tf)={0};
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx'],data, 'Sheet1', 'A1');
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx'],Header, 'Sheet1', 'A1');
                clear data;
                cd ..\..\..
                cd('Results_matched');
                loc= location(ColumnIndexMergeFile);
                out=regexp(loc,'\','split');
                out=out{1};
                length(out);
                CurrentPeptide=out{length(out)-2};
                CurrentReplicate=out{length(out)-1};
                cd(CurrentPeptide);
                cd(CurrentReplicate);
                directories{idxx}=CurrentPeptide;
                idxx=idxx+1;
                delete([strcat(CurrentPeptide, CurrentReplicate),'.xlsx']);
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],FileMat2, 'Sheet1', 'A1');
                Header = string(FileMat2(1,:));
                Header = cell(1,length(Header));
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],Header, 'Sheet1', 'A1');
                cd ..\..\..

            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%% replicate 2
for RowIndexMergeFile= 2:3: length(File)
    for ColumnIndexMergeFile= 2:3:length(File)
        [num, txt, FileMat] = xlsread(string(location(RowIndexMergeFile)));
        if RowIndexMergeFile== ColumnIndexMergeFile
            %do nothing
        elseif ColumnIndexMergeFile>RowIndexMergeFile
            [num, txt, FileMat2] = xlsread(string(location(ColumnIndexMergeFile)));
            %matching files
            logicalArray=ismember(string(FileMat2(1,:)),string(FileMat(1,:)));
            num=0;
            num2=0;
            % if Header is same
            if    any (logicalArray(:),1)
                % find  the index of blocks of 2 files
                [row_ind1, row_end1] = XRayDosageDataBlockIndicescaseOne(FileMat);
                [row_ind2, row_end2] = XRayDosageDataBlockIndicesCaseTwo(FileMat2);

                for pos=1:length(row_ind1)
                    %Extracting rows between the index of first row and the last row within a block of X-ray dosage
                    row= row_ind1(pos): row_end1(pos);
                    length(row);
                    %%%%%%%%%%%%%%%%%%%%%%%;
                    row2= row_ind2(pos): row_end2(pos);
                    index=5;
                    data(row_ind1(pos)+num: row_end1(pos)+num,1:4)= FileMat(row_ind1(pos): row_end1(pos),1:4);
                    data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,1:4)= FileMat2(row_ind2(pos): row_end2(pos),1:4);
                    FileMat2(row_ind2(pos): row_end2(pos),1:4)= {0};
                    % add the length
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% for residues
                    for index2= 5: size(FileMat2,2)

                        if string(FileMat2(1,index2))== string(FileMat(1,index))
                            data(row_ind1(pos)+num: row_end1(pos)+num,index2)= FileMat(row_ind1(pos): row_end1(pos),index);
                            data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,index2)= FileMat2(row_ind2(pos): row_end2(pos),index2);
                            FileMat2(row_ind2(pos): row_end2(pos),index2)={0};
                            index=min(size(FileMat,2),index+1);
                        else
                            data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,index2)= FileMat2(row_ind2(pos): row_end2(pos),index2);
                            FileMat2(row_ind2(pos): row_end2(pos),index2)={0};
                        end
                    end

                    num=num+length(row2);
                    num2=num2+length(row);
                end
                FirstColumn = data(:,1);
                Storage = [];

                for index = 1:size(FirstColumn,1)
                    CellTemp = FirstColumn(index);
                    if iscellstr(CellTemp)
                        strTemp = string(CellTemp);
                        Storage = [Storage; [strTemp index] ];
                    end
                end

                idx=[];
                jidx=1;
                for id= 2:2:size(Storage,1)
                    idx(jidx)= Storage(id,2);
                    jidx=jidx+1;
                end
                data(idx,:)=[];

                if ~isfolder(strcat(pwd,'\\Resultsnew1'))
                    mkdir('Resultsnew1');
                end

                cd ('Resultsnew1');

                Header = string(FileMat2(1,:));
                loc= location(RowIndexMergeFile);
                out=regexp(loc,'\','split');
                out=out{1};
                length(out);
                CurrentPeptide=out{length(out)-2};
                CurrentReplicate=out{length(out)-1};
                mkdir(CurrentPeptide);
                cd(CurrentPeptide);
                mkdir(CurrentReplicate);
                cd(CurrentReplicate);
                delete([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx']);
                temp = data(:,5: size(data,2));
                temp(cellfun(@(x)isscalar(x) && isnan(x),temp)) = {0};
                data(:,5: size(data,2)) = temp;
                tf=cellfun('isempty',data );
                data(tf)={0};
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx'],data, 'Sheet1', 'A1');
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx'],Header, 'Sheet1', 'A1');
                clear data;
                cd ..\..\..
                cd('Results_matched');
                loc= location(ColumnIndexMergeFile);
                out=regexp(loc,'\','split');
                out=out{1};
                length(out);
                CurrentPeptide=out{length(out)-2};
                CurrentReplicate=out{length(out)-1};
                cd(CurrentPeptide);
                cd(CurrentReplicate);
                delete([strcat(CurrentPeptide, CurrentReplicate),'.xlsx']);
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],FileMat2, 'Sheet1', 'A1');
                Header = string(FileMat2(1,:));
                Header = cell(1,length(Header));
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],Header, 'Sheet1', 'A1');
                cd ..\..\..

            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%% replicate 3
for RowIndexMergeFile= 3:3: length(File)
    for ColumnIndexMergeFile= 3:3:length(File)
        [num, txt, FileMat] = xlsread(string(location(RowIndexMergeFile)));
        if RowIndexMergeFile== ColumnIndexMergeFile
            %do nothing

        elseif ColumnIndexMergeFile>RowIndexMergeFile
            [num, txt, FileMat2] = xlsread(string(location(ColumnIndexMergeFile)));
            %matching files
            logicalArray=ismember(string(FileMat2(1,:)),string(FileMat(1,:)));
            num=0;
            num2=0;
            % if Header is same
            if    any (logicalArray(:),1)
                % find  the index of blocks of 2 files
                [row_ind1, row_end1] = XRayDosageDataBlockIndicesCaseOne(FileMat);
                [row_ind2, row_end2] = XRayDosageDataBlockIndicesCaseTwo(FileMat2);

                for pos=1:length(row_ind1)
                    %Extracting rows between the index of first row and the last row within a block of X-ray dosage
                    row= row_ind1(pos): row_end1(pos);
                    length(row);
                    %%%%%%%%%%%%%%%%%%%%%%%
                    row2= row_ind2(pos): row_end2(pos);
                    index=5;
                    data(row_ind1(pos)+num: row_end1(pos)+num,1:4)= FileMat(row_ind1(pos): row_end1(pos),1:4);
                    data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,1:4)= FileMat2(row_ind2(pos): row_end2(pos),1:4);
                    FileMat2(row_ind2(pos): row_end2(pos),1:4)= {0};
                    % add the length
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% for residues

                    for index2= 5: size(FileMat2,2)


                        if string(FileMat2(1,index2))== string(FileMat(1,index))
                            data(row_ind1(pos)+num: row_end1(pos)+num,index2)= FileMat(row_ind1(pos): row_end1(pos),index);
                            data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,index2)= FileMat2(row_ind2(pos): row_end2(pos),index2);
                            FileMat2(row_ind2(pos): row_end2(pos),index2)={0};
                            index=min(size(FileMat,2),index+1);

                        else

                            data(row_ind2(pos)+length(row)+num2: row_end2(pos)+length(row)+num2,index2)= FileMat2(row_ind2(pos): row_end2(pos),index2);
                            FileMat2(row_ind2(pos): row_end2(pos),index2)={0};
                        end
                    end
                    num=num+length(row2);
                    num2=num2+length(row);

                end
                FirstColumn = data(:,1);
                Storage = [];

                for index = 1:size(FirstColumn,1)
                    CellTemp = FirstColumn(index);
                    if iscellstr(CellTemp)
                        strTemp = string(CellTemp);
                        Storage = [Storage; [strTemp index] ];
                    end
                end
                idx=[];
                jidx=1;
                for id= 2:2:size(Storage,1)
                    idx(jidx)= Storage(id,2);
                    jidx=jidx+1;
                end
                data(idx,:)=[];
                if ~isfolder(strcat(pwd,'\\Resultsnew1'))
                    mkdir('Resultsnew1');
                end

                cd ('Resultsnew1');

                Header = string(FileMat2(1,:));
                loc= location(RowIndexMergeFile);
                out=regexp(loc,'\','split');
                out=out{1};
                length(out);
                CurrentPeptide=out{length(out)-2};
                CurrentReplicate=out{length(out)-1};
                mkdir(CurrentPeptide);
                cd(CurrentPeptide);
                mkdir(CurrentReplicate);
                cd(CurrentReplicate);
                delete([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx']);
                temp = data(:,5: size(data,2));
                temp(cellfun(@(x)isscalar(x) && isnan(x),temp)) = {0};
                data(:,5: size(data,2)) = temp;
                tf=cellfun('isempty',data );
                data(tf)={0};
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx'],data, 'Sheet1', 'A1');
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'o.xlsx'],Header, 'Sheet1', 'A1');
                clear data;
                cd ..\..\..
                cd('Results_matched');
                loc= location(ColumnIndexMergeFile);
                out=regexp(loc,'\','split');
                out=out{1};
                length(out);
                CurrentPeptide=out{length(out)-2};
                CurrentReplicate=out{length(out)-1};
                cd(CurrentPeptide);
                cd(CurrentReplicate);
                delete([strcat(CurrentPeptide, CurrentReplicate),'.xlsx']);
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],FileMat2, 'Sheet1', 'A1');
                Header = string(FileMat2(1,:));
                Header = cell(1,length(Header));
                xlswrite([strcat(CurrentPeptide, CurrentReplicate),'.xlsx'],Header, 'Sheet1', 'A1');
                cd ..\..\..

            end
        end
    end
end
if isfolder(strcat(pwd,'\\Resultsnew1'))
    cd('Results_matched');
    %% removing the merge directories

    for RowIndexMergeFile=1:size(directories,2)
        rmdir(string(directories(RowIndexMergeFile)),'s');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % saving the unmodified folders in new Directories
    cd ..
    direc= char(strcat(pwd,strcat('\Results_matched')));
    copyfile(direc,'Resultsnew1');
else
    % do nothing
end




