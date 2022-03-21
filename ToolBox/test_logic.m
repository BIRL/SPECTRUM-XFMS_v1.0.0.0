[num, txt, File] = xlsread('datamock.xlsx')
File = string(File)
File=File(:,any(~ismissing(File(2:length(File),:))))
% Step 3: Extract indices of data blocks at various X Ray dosages
[row_ind, row_end] = XRayDosageDataBlockIndices(File);

% Step 4: Extract the  names of oxidized residues
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
    total_1 = sum(File_Main(row,2));
    % file_main(row,15) contains  number of rows in column 3 - oxidized area(R1)
    total_2 = sum(File_Main(row,3));
    % sum of columns - oxidized area(R1) and Unoxidized area(R1)
    R1_total=total_1+ total_2;
    % sum is saved in 1st column of 'Result1' matrix
    Result1(pos,1)= R1_total;



    % Step 7: Calculating F values for each oxidized residue in a peptide
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%% for Replicate 1  %%%%%%%%%%%%%%%%%%%%%%
    File_Main2=zeros(length(File_Main),16)
    File_Main2(:,Index_var)=File_Main(:,Index_var)
    % Defining area as column 3 - oxidized area(R1)
    Area= File_Main(:,3);
    % Transpose of variable 'area' is taken for computataion purposes
    Area=Area';

    % column 22 contains information about overlapping of retention time of
    % 2 residues
    col= File_Main(:,4);
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
        for i=row(pos): row_end(pos)
            % checking wheter the residue has overlap retention time with
            % any other residues
            if col(i)==1 && File_Main2(i,j)== 1
                %Extracting the column name to Get the name of the amino acid
                %and its reacticity.
                name_col= convertStringsToChars(colname(j));
                Amino_acid1=name_col(1);
                IndexofResidue = strfind(amino_acids,Amino_acid1);
                RC_1a= RC_1a + reactivity(IndexofResidue);
                % checking the neighbouring columns to identify the overlapped residue
                if File_Main2(i,j+1)==1
                    %Extracting the column name to Get the name of the amino acid
                    %and its reacticity.
                    name_col= convertStringsToChars(colname(j+1));
                    Amino_acid2=name_col(1);
                    IndexofResidue = strfind(amino_acids,Amino_acid2);
                    RC_1b= RC_1b + reactivity(IndexofResidue);
                    % checking the neighbouring columns to identify the overlapped residue
                else if File_Main2(i,j-1)==1
                        %Extracting the column name to Get the name of the amino acid
                        %and its reacticity.
                        name_col= convertStringsToChars(colname(j-1));
                        Amino_acid2=name_col(1);
                        IndexofResidue = strfind(amino_acids,Amino_acid2);
                        RC_1b= RC_1b + reactivity(IndexofResidue);

                        % checking the neighbouring columns to identify the overlapped residue
                else  if File_Main2(i,j+2)==1
                        %Extracting the column name to Get the name of the amino acid
                        %and its reacticity.
                        name_col= convertStringsToChars(colname(j+2));
                        Amino_acid2=name_col(1);
                        IndexofResidue = strfind(amino_acids,Amino_acid2);
                        RC_1b= RC_1b + reactivity(IndexofResidue);
                        % checking the neighbouring columns to identify the overlapped residue
                else if File_Main2(i,j-2)==1
                        %Extracting the column name to Get the name of the amino acid
                        %and its reacticity.
                        name_col= convertStringsToChars(colname(j-2));
                        Amino_acid2=name_col(1);
                        IndexofResidue = strfind(amino_acids,Amino_acid2);
                        RC_1b= RC_1b + reactivity(IndexofResidue);
                else if File_Main2(i,j+3)==1
                        %Extracting the column name to Get the name of the amino acid
                        %and its reacticity.
                        name_col= convertStringsToChars(colname(j+3));
                        Amino_acid2=name_col(1);
                        IndexofResidue = strfind(amino_acids,Amino_acid2);
                        RC_1b= RC_1b + reactivity(IndexofResidue);
                else if File_Main2(i,j-3)==1
                        %Extracting the column name to Get the name of the amino acid
                        %and its reacticity.
                        name_col= convertStringsToChars(colname(j-3));
                        Amino_acid2=name_col(1);
                        IndexofResidue = strfind(amino_acids,Amino_acid2);
                        RC_1b= RC_1b + reactivity(IndexofResidue);
                end
                end
                end
                end
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
            else if col(i)==2 && File_Main2(i,j)== 1
                    %Extracting the column name to Get the name of the amino acid
                    %and its reacticity.
                    name_col= convertStringsToChars(colname(j));
                    Amino_acid1=name_col(1);
                    IndexofResidue = strfind(amino_acids,Amino_acid1);
                    RC_1a= RC_1a + reactivity(IndexofResidue);
                    % checking the neighbouring columns to identify the overlapped residue
                    if File_Main2(i,j+1)==2
                        name_col1= convertStringsToChars(colname(j+1));
                        Amino_acid2a=name_col1(1);
                        IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                        RC_1b= RC_1b + reactivity(IndexofResidue1);
                        if File_Main2(i,j-1)==1
                            name_col2= convertStringsToChars(colname(j-1));
                            Amino_acid2b=name_col2(1);
                            IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                            RC_1c= RC_1c + reactivity(IndexofResidue2);
                        else if File_Main2(i,j-2)==1
                                name_col2= convertStringsToChars(colname(j-2));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                        else if File_Main2(i,j-3)==1
                                name_col2= convertStringsToChars(colname(j-3));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                        else if File_Main2(i,j+2)==1
                                name_col2= convertStringsToChars(colname(j+2));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                        else if File_Main2(i,j+3)==1
                                name_col2= convertStringsToChars(colname(j+3));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                        end
                        end
                        end
                        end
                        end

                    else if File_Main2(i,j-1)==2
                            name_col1= convertStringsToChars(colname(j-1));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-2)==1
                                    name_col2= convertStringsToChars(colname(j-2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-3)==1
                                    name_col2= convertStringsToChars(colname(j-3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+2)==1
                                    name_col2= convertStringsToChars(colname(j+2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+3)==1
                                    name_col2= convertStringsToChars(colname(j+3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end

                    else if File_Main2(i,j-2)==2
                            name_col1= convertStringsToChars(colname(j-2));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-1)==1
                                    name_col2= convertStringsToChars(colname(j-1));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-3)==1
                                    name_col2= convertStringsToChars(colname(j-3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+2)==1
                                    name_col2= convertStringsToChars(colname(j+2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+3)==1
                                    name_col2= convertStringsToChars(colname(j+3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end
                    else if File_Main2(i,j+2)==2
                            name_col1= convertStringsToChars(colname(j+2));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-1)==1
                                    name_col2= convertStringsToChars(colname(j-1));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-3)==1
                                    name_col2= convertStringsToChars(colname(j-3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-2)==1
                                    name_col2= convertStringsToChars(colname(j-2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+3)==1
                                    name_col2= convertStringsToChars(colname(j+3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end
                    else if File_Main2(i,j+3)==2
                            name_col1= convertStringsToChars(colname(j+3));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-1)==1
                                    name_col2= convertStringsToChars(colname(j-1));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-3)==1
                                    name_col2= convertStringsToChars(colname(j-3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-2)==1
                                    name_col2= convertStringsToChars(colname(j-2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+2)==1
                                    name_col2= convertStringsToChars(colname(j+2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end

                    else if File_Main2(i,j-3)==2
                            name_col1= convertStringsToChars(colname(j-3));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-1)==1
                                    name_col2= convertStringsToChars(colname(j-1));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+3)==1
                                    name_col2= convertStringsToChars(colname(j+3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-2)==1
                                    name_col2= convertStringsToChars(colname(j-2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+2)==1
                                    name_col2= convertStringsToChars(colname(j+2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end


                    end
                    end
                    end
                    end
                    end
                    end
                    % sum of reactivitity constants of 2 aminoacids
                    total_RC1= RC_1b+ RC_1a;
                    total_RC2= RC_1b+ RC_1c;
                    RC=total_RC1/(total_RC1+total_RC2);
                    fprintf('Reactivity of %s with %s and %s at row %f is %f .\n ',Amino_acid1,Amino_acid2a,Amino_acid2b, i,RC)
                    % Divided area of residues having same Retention time
                    Total_area_R=Area(i)* RC;
                    % Add the area to the variable count
                    count=count+Total_area_R;

            else if col(i)==2 && File_Main2(i,j)== 2
                    %Extracting the column name to Get the name of the amino acid
                    %and its reacticity.
                    name_col= convertStringsToChars(colname(j));
                    Amino_acid1=name_col(1);
                    IndexofResidue = strfind(amino_acids,Amino_acid1);
                    RC_1a= RC_1a + reactivity(IndexofResidue);
                    if File_Main2(i,j+1)==1
                        name_col1= convertStringsToChars(colname(j+1));
                        Amino_acid2a=name_col1(1);
                        IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                        RC_1b= RC_1b + reactivity(IndexofResidue1);
                        if File_Main2(i,j-1)==1 || File_Main2(i,j-1)==2
                            name_col2= convertStringsToChars(colname(j-1));
                            Amino_acid2b=name_col2(1);
                            IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                            RC_1c= RC_1c + reactivity(IndexofResidue2);
                        else if File_Main2(i,j-2)==1 || File_Main2(i,j-2)==2
                                name_col2= convertStringsToChars(colname(j-2));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                        else if File_Main2(i,j-3)==1 || File_Main2(i,j-3)==2
                                name_col2= convertStringsToChars(colname(j-3));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                        else if File_Main2(i,j+2)==1 || File_Main2(i,j+2)==2
                                name_col2= convertStringsToChars(colname(j+2));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                        else if File_Main2(i,j+3)==1 || File_Main2(i,j+3)==2
                                name_col2= convertStringsToChars(colname(j+3));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                        end
                        end
                        end
                        end
                        end

                    else if File_Main2(i,j-1)==1
                            name_col1= convertStringsToChars(colname(j-1));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1 || File_Main2(i,j+1)==2
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-2)==1 || File_Main2(i,j-2)==2
                                    name_col2= convertStringsToChars(colname(j-2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-3)==1 || File_Main2(i,j-3)==2
                                    name_col2= convertStringsToChars(colname(j-3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+2)==1 || File_Main2(i,j+2)==2
                                    name_col2= convertStringsToChars(colname(j+2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+3)==1 || File_Main2(i,j+3)==2
                                    name_col2= convertStringsToChars(colname(j+3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end

                    else if File_Main2(i,j-2)==1
                            name_col1= convertStringsToChars(colname(j-2));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1 || File_Main2(i,j+1)==2
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-1)==1 || File_Main2(i,j-1)==2
                                    name_col2= convertStringsToChars(colname(j-1));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-3)==1 || File_Main2(i,j-3)==2
                                    name_col2= convertStringsToChars(colname(j-3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+2)==1 || File_Main2(i,j+2)==2
                                    name_col2= convertStringsToChars(colname(j+2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+3)==1 || File_Main2(i,j+3)==2
                                    name_col2= convertStringsToChars(colname(j+3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end
                    else if File_Main2(i,j+2)==1
                            name_col1= convertStringsToChars(colname(j+2));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1 || File_Main2(i,j+1)==2
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-1)==1 || File_Main2(i,j-1)==2
                                    name_col2= convertStringsToChars(colname(j-1));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-3)==1 || File_Main2(i,j-3)==2
                                    name_col2= convertStringsToChars(colname(j-3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-2)==1 || File_Main2(i,j-2)==2
                                    name_col2= convertStringsToChars(colname(j-2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+3)==1 || File_Main2(i,j+3)==2
                                    name_col2= convertStringsToChars(colname(j+3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end
                    else if File_Main2(i,j+3)==1
                            name_col1= convertStringsToChars(colname(j+3));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1 || File_Main2(i,j+1)==2
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-1)==1 || File_Main2(i,j-1)==2
                                    name_col2= convertStringsToChars(colname(j-1));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-3)==1 || File_Main2(i,j-3)==2
                                    name_col2= convertStringsToChars(colname(j-3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-2)==1 || File_Main2(i,j-2)==2
                                    name_col2= convertStringsToChars(colname(j-2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+2)==1 || File_Main2(i,j+2)==2
                                    name_col2= convertStringsToChars(colname(j+2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end

                    else if File_Main2(i,j-3)==1
                            name_col1= convertStringsToChars(colname(j-3));
                            Amino_acid2a=name_col1(1);
                            IndexofResidue1 = strfind(amino_acids,Amino_acid2a);
                            RC_1b= RC_1b + reactivity(IndexofResidue1);
                            if File_Main2(i,j+1)==1  || File_Main2(i,j+1)==2
                                name_col2= convertStringsToChars(colname(j+1));
                                Amino_acid2b=name_col2(1);
                                IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-1)==1 || File_Main2(i,j-1)==2
                                    name_col2= convertStringsToChars(colname(j-1));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+3)==1 || File_Main2(i,j+1)==2
                                    name_col2= convertStringsToChars(colname(j+3));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j-2)==1 || File_Main2(i,j-2)==2
                                    name_col2= convertStringsToChars(colname(j-2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            else if File_Main2(i,j+2)==1 || File_Main2(i,j+2)==2
                                    name_col2= convertStringsToChars(colname(j+2));
                                    Amino_acid2b=name_col2(1);
                                    IndexofResidue2 = strfind(amino_acids,Amino_acid2b);
                                    RC_1c= RC_1c + reactivity(IndexofResidue2);
                            end
                            end
                            end
                            end
                            end


                    end
                    end
                    end
                    end
                    end
                    end
                    % sum of reactivitity constants of 2 aminoacids
                    total_RC1a= RC_1a+ RC_1b;
                    total_RC2a= RC_1a+ RC_1c;
                    RCa=total_RC1/(total_RC1+total_RC2);
                    RCb=total_RC2/(total_RC1+total_RC2);
                    RC=RCa+RCb
                    fprintf('Reactivity of %s with %s and %s at row %f is %f .\n ',Amino_acid1,Amino_acid2a,Amino_acid2b, i,RC)
                    % Divided area of residues having same Retention time
                    Total_area_R=Area(i)* RC;
                    % Add the area to the variable count
                    count=count+Total_area_R;
                    % Sum of area of  residue that doesnot operlap with other residue






                    % Sum of area of  residue that doesnot operlap with other residue

            else if  col(i)==0 && File_Main2(i,j)== 1
                    count=count+Area(i);
            end
            end
            end
            end
            RC_1a=0;
            RC_1b=0;
            RC_1c=0;
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
