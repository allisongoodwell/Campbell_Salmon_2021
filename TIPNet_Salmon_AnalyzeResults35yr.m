%%%%%% TIPNet for CRB Salmon study: Results analysis %%%%%
% Allison Goodwell and Nick Campbell
% June 2020
% This code loads in results files from TIPNet_Salmon_35yr.m
% and makes some plots of results

% notes: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

load('RESULTS_35yr.mat')

%15 columns, column 1 is the date
data = xlsread('masterfile35yr.xlsx');

year = data(:,1);
salmon_pops = data(:,2:7);
flows = data(:,8:12);
watertempdays = data(:,13:18); 
PDO = data(:,19);
PPT = data(:,20);
airtemp = data(:,21);
degdays = data(:,22);
SWE = data(:,23:26);

salmon_pops_orig = salmon_pops;

%AEG edit 10/5/20: keep BON, but for others - fraction of BON
for i =2:6

    vect=salmon_pops(:,i)./ salmon_pops(:,1); %
    salmon_pops_orig(:,i)=vect;
    %omit outliers
    minval = mean(vect)-1.5*iqr(vect);
    maxval = mean(vect)+1.5*iqr(vect);
    vect(vect>maxval)=maxval;
    vect(vect<minval)=minval;
    salmon_pops(:,i)=vect;
    
end

%%

%AEG edit 10/15/20: set outliers in driving variables to max values
for i = 1:size(watertempdays,2)
    vect = watertempdays(:,i);
    maxval = mean(vect)+1.5*iqr(vect);
    vect(vect>maxval)=maxval;
    watertempdays(:,i)=vect;   
end

maxPPT = mean(PPT)+1.5*iqr(PPT);
PPT(PPT>maxPPT)=maxPPT;

minairtemp = mean(airtemp)-1.5*iqr(airtemp);
airtemp(airtemp<minairtemp)=minairtemp;


AllData = [salmon_pops, flows, watertempdays, PDO, PPT, airtemp, degdays, SWE];

fishsitenames = {'BON','MCN','IHR','LWG','PRD','WEL'};
flownames = {'QBON','Q_LWG','QPRD','QDWR', 'SpillDWR'};
tempnames = {'TempDaysJDA68','TempDaysJDA70','TempDaysLWG68',...
    'TempDaysLWG70','TempDaysPRD68','TempDaysPRD70'};
PDO_name = {'PDO'};
PPT_name = {'Precip'};
Airtemp_name = {'Ta'};
DegDay_name = {'TaDays90'};
SWE_name = {'SWESnake','SWEMidCRB','SWEUpCRB','SWELowCRB'};

varnames_nofish = [flownames tempnames PDO_name PPT_name Airtemp_name DegDay_name, SWE_name];

%plot some data
figure(1)
subplot(6,1,1)
plot(salmon_pops(:,2:6))
legend(fishsitenames)
subplot(6,1,2)
plot(flows)
legend(flownames)
subplot(6,1,3)
plot(watertempdays)
legend(tempnames)
subplot(6,1,4)
plot(PDO)
legend('PDO')
subplot(6,1,5)
plot(SWE)
legend(SWE_name)
subplot(6,1,6)
plot(PPT)
legend('precip')


Ilags = ent.I_tau_normbyH; %lags, source, target

markertypes = ['.','.','.','.','.','.',... 
    '*','*','*','*','*',...
    'o','o','o','o','o','o',...
    'p','s','d','d','d','d','d','d'];
colorvect = [jet(6);jet(5);jet(7);jet(7)];

figure(2)
for i =1:6
    subplot(3,2,i)
    hold on
    for s = 7:size(AllData,2)
       plot(0:7,Ilags(:,s,i),'Marker',markertypes(s),'MarkerSize',15,'LineStyle','none','Color',colorvect(s,:),'MarkerFaceColor',colorvect(s,:)); 
        [maxval maxind]=max(Ilags(:,s,i));
    
        MaxImatrix(i,s-6)=maxval;
        MaxIlag(i,s-6)=maxind-1;
        ax=gca;
        ax.FontSize = 13;
    
    end
    
    xlim([0 7])
    ylim([.1 1])
    xlabel('time lag (years)')
    ylabel('normalized MI')
    if i==1
    %legend(varnames_nofish,'Location','EastOutside')
    end
    title(fishsitenames(i))
    
end





%% do a lagged correlation analysis (pasting code from Alyssa's MS report...)

drivers_indices = 1:size(AllData,2);
targets_indices = 1:6;
nlags = 7;
 
[Corrvals, pvals] = corrcoef(AllData,'Rows','pairwise');
  
for t = 0:nlags 
 
drivers_lagged = AllData(1:end-t, drivers_indices); 
targets_current = AllData((t+1):end, targets_indices); 
new_variables = [drivers_lagged targets_current]; 

[Corrvals_lagged(:,:,t+1), pvals_lagged(:,:,t+1)] = corrcoef(new_variables,'Rows','pairwise');
 
end

k=1;


for  j = (size(AllData,2)+1):(size(AllData,2)+6) %targets (salmon)
    
    for i = 7:size(AllData,2) %lagged sources (all variables)
        
        vector=reshape(Corrvals_lagged(j,i,1:nlags+1),[1,nlags+1]);
        vectorp=reshape(pvals_lagged(j,i,1:nlags+1),[1,nlags+1]);
        vector(vectorp > 0.005) = nan;
        numnans = length(find(isnan(vector)));
        
        [maxval maxind]=max(vector);
    
        MaxCorrmatrix(j-size(AllData,2),i-6)=maxval;
        MaxCorrlag(j-size(AllData,2),i-6)=maxind-1;
    
        figure(3)
        subplot(3,2,k)
        hold on
        plot(0:7, vector,'Marker',markertypes(i),'MarkerSize',15,'Color',colorvect(i,:),'LineStyle','none','MarkerFaceColor',colorvect(i,:));
        %yline(0,'--r');
        title (sprintf('%s',varnames_all{j-size(AllData,2)}));
        xlabel('time lag (years)')
        %sgtitle('Lagged Correlations');
        ax = gca;
        ax.FontSize = 13;
        ylim([-1 1])
        xlim([0 7])
        
        

    end
    
    k=k+1;
end


%% make some plots of cases that find lagged MI but not correlation

%TempDays PRD (68 and 70) to IHR and WEL at time lags of 6-7 years

BON = salmon_pops(:,1);
IHR_frac =salmon_pops(:,3);
WEL_frac = salmon_pops(:,6);
PRD_frac = salmon_pops(:,5);

BON = (BON-min(BON))./range(BON);
WEL_frac = (WEL_frac-min(WEL_frac))./range(WEL_frac);
PRD_frac = (PRD_frac-min(PRD_frac))./range(PRD_frac);

SWE_midCRB = SWE(:,2);
SWE_midCRB = (SWE_midCRB-min(SWE_midCRB))./range(SWE_midCRB);

PPT = (PPT-min(PPT))./range(PPT);
PDO = (PDO-min(PDO))./range(PDO);

TD_LWG68=watertempdays(:,3);
TD_LWG68 = (TD_LWG68-min(TD_LWG68))./range(TD_LWG68);

Spill = flows(:,5);
Spill = (Spill-min(Spill))./range(Spill);


figure(4) %lagged sources to BON
subplot(2,3,1)
plot(PDO(1:end-3),BON(4:end),'.b','MarkerSize',10)
title('PDO to BON')
xlabel('PDO(t-3)')
ylabel('BON(t)')

subplot(2,3,2)
plot(SWE_midCRB(1:end-4),BON(5:end),'.b','MarkerSize',10)
title('SWE to BON')
xlabel('SWE(t-4)')
ylabel('BON(t)')
subplot(2,3,3)
plot(PPT(1:end-6),BON(7:end),'.b','MarkerSize',10)
title('Precip to BON')
xlabel('Precip(t-6)')
ylabel('BON(t)')
 %targets of TempDays LWG 68
subplot(2,3,4)
plot(TD_LWG68(1:end),WEL_frac(1:end),'.b','MarkerSize',10)
title('TDLWG68 to WEL/BON')
ylabel('WEL(t)/BON(t)')
xlabel('TDLWG68(t)')
subplot(2,3,5)
plot(TD_LWG68(1:end-4),BON(5:end),'.b','MarkerSize',10)
title('TDLWG68 to BON')
xlabel('TDLWG68(t-4)')
ylabel('BON(t)')
subplot(2,3,6)
plot(Spill(1:end),PRD_frac(1:end),'.b','MarkerSize',10)
title('Spill to PRD/BON')
xlabel('Spill(t)')
ylabel('PRD(t)/BON(t)')

%% Circos plots for sources...correlation vs MI

ntargets=6;

foldername = 'CircosFiles';
if exist(foldername)~=7
    mkdir(foldername);
end
addpath(foldername);

colR = round(colorvect.*255);

for j =1:2

    if j==1
        fileID = fopen([foldername,'/MI35yr','.txt'],'w+');
        matrix = MaxImatrix;
        colSource = 0;
    elseif j==2
        fileID = fopen([foldername,'/Corr35yr','.txt'],'w+');
        matrix = abs(MaxCorrmatrix);
        matrix(isnan(matrix))=0;
        colSource = 255;
    end
    
    %print first line: label label 1 2 3 4.....nvars
    fprintf(fileID,'label label ');
    for v =1:ntargets
        fprintf(fileID,'%d ',v);
    end
    fprintf(fileID,'\n');
    
    %print second line: label label colors (all same for source vars)
    fprintf(fileID,'label label ');
    for v =1:ntargets
        fprintf(fileID,'%d,%d,%d ',colSource,colSource,colSource);
    end
    fprintf(fileID,'\n'); 
    
    %print third line: label label target varnames
    fprintf(fileID,'label label ');
    for v =1:ntargets
        fprintf(fileID,'%s ',fishsitenames{v});
    end
    fprintf(fileID,'\n');
    
    %print 3rd -- lines: color sourcenames
    
    for vp = 1:size(matrix,2)
        fprintf(fileID,'%d,%d,%d %s ',colR(vp,1), colR(vp,2), colR(vp,3), varnames_nofish{vp});

        for v =1:ntargets
            fprintf(fileID,'%d  ',round(matrix(v,vp).*100));
        end

        fprintf(fileID,'\n');
    end


end

