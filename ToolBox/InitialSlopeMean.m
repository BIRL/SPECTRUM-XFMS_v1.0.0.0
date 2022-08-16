%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM-XFMS: A MATLAB Toolbox to Analyze X-ray Footprinting Mass Spectrometry Data %
%                                    Version 1.0.0                                     %
%        Copyright (c) Biomedical Informatics & Engineering Research Laboratory,       %
%           Lahore University of Management Sciences Lahore (LUMS), Pakistan.          %
%                           (http://biolabs.lums.edu.pk/BIRL)                          %
%                                (safee.ullah@gmail.com)                               %
%                            Last Modified on: 26-July-2022                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function reads the  FileMean generated in the Step 11.
% This file includes the mean and standard error of the residues in three repliactes
% It fit the exponential curve on the data of each residues to creat a plot
function [ExponentialFit, DataSlope] =  InitialSlopeMean( FullFileName3,header,FileMean,File,col_names)
Dosage= [0;10;25;50;75;100];
Idx_row= 1;
Index_col=1;
Index_row=1;
   IdxGradient=1;
    TotalColumn= size(FileMean);
    TotalColumn=TotalColumn(2);
% for each column of (1-F) in result
for j= 2:3:TotalColumn
    % For each number of rows
    for i=1:length(Dosage)
        %Itratively Extracting the (1-F) and dosage for a single residue
        Data_Mean(i,:)=[FileMean(i,j), Dosage(i,1),FileMean(i,j+1),FileMean(i,j+2)];
    end
    % converting it to matrix for easy plotting
  Data_Mean=cell2mat(Data_Mean);
    % Defining the Y and X variables as Dosage and '1-F' respectively
    x= Data_Mean(:,2);
    y= Data_Mean(:,1);
    error_1=Data_Mean(:,3);
    weight=Data_Mean(:,4);
   
    % fitting the exponential curve
    %%%% Loop to remove Zero value from fit at higer dosage
    FilteredY=y;
    FilteredX=x;
    FilteredWeight= weight;
    for IndexRemove =2: length(y)
           if FilteredY(IndexRemove) == 1
              FilteredY(IndexRemove) = nan;
                FilteredX(IndexRemove) = nan;
                FilteredWeight(IndexRemove) = nan;
           end
    end
    FilteredY = rmmissing(FilteredY);
    FilteredX = rmmissing(FilteredX);
    FilteredWeight = rmmissing(FilteredWeight);
    ExponentialFit=fit(x,y,'exp1');
    WeightedExpFit=fit(FilteredX,FilteredY,ExponentialFit,'Weight',FilteredWeight);
%%%%%%intervals
if length(FilteredY)< 3
figure (1);
    set(gcf,'Visible', 'off');
    % plotting the curve
    % plotting confidence interval in graph

plot(WeightedExpFit);
   errorbar(x ,y,error_1,'ko','both');
hold on 
plot(WeightedExpFit );
plot(x,y,'.','MarkerSize',15);
hold off

else

int=confint(WeightedExpFit);
int(3,:)=coeffvalues(WeightedExpFit);

    figure (1);
    set(gcf,'Visible', 'off');
    % plotting the curve
    % plotting confidence interval in graph
y_fit = int(3,1)*exp(int(3,2)*x);
y_cbl = int(1,1)*exp(int(1,2)*x);
y_cbh = int(2,1)*exp(int(2,2)*x);
plot(x,y_fit,'r',x,y_cbl,'b:',x,y_cbh,'b:');
   errorbar(x ,y,error_1,'ko','both');
hold on 
plot(WeightedExpFit );
plot(x,y,'.','MarkerSize',15);

plot(x,y_cbl,'b:',x,y_cbh,'b:');

hold off
end 
    % Naming the plot as the name of the residue

    name_x= col_names(IdxGradient);
 cd ('Result');
    cd ( FullFileName3);
    title_string = strcat(name_x, ' -Mean Exponential Fit Curve');
    title(title_string, 'Fontsize', 20 , 'fontweight', 'bold', 'Color', [0.8 0 0]);
    xlabel('Dose', 'Fontsize', 16, 'fontweight', 'bold', 'Color', [0 0 0]);
    ylabel('1 - F', 'Fontsize', 16, 'fontweight', 'bold', 'Color', [0 0 0]);
      legend({'Errorbar','Fit','Data', 'Interval'},'Location','southwest');
       xlim([-5 105])
    plotname = strcat(name_x,'_fit_Mean.png');
    saveas(gcf,plotname);
    % Getting the coefficient values

    coeff=coeffvalues(ExponentialFit);
    slope=coeff(2);
    % saving the slope Value in a matrix Data_slope
    DataSlope(Idx_row: Idx_row+2 ,Index_col:Index_col+1)= int;
     Gra= int(3,1)*int(3,2)*exp(int(3,2)*y(1));
    DataSlope(Idx_row+4 ,Index_col+1)= Gra;
    Index_col=Index_col+2;
    Index_row=Index_row+1;
    IdxGradient=IdxGradient+1;
    clear Data_Mean
    cd ..\..
end