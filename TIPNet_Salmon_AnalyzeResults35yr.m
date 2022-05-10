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
plot(TD_LWG68(1:end-3),BON(4:end),'.b','MarkerSize',10)
title('TDLWG68 to BON')
xlabel('TDLWG68(t-3)')
ylabel('BON(t)')

subplot(2,3,6)
plot(PDO(1:end-5),WEL_frac(6:end),'.b','MarkerSize',10)
title('PDO to WEL')
xlabel('PDO(t-5)')
ylabel('BON(t)')



%% do a lagged correlation analysis for comparison with MI

drivers_indices = 1:size(AllData,2);
targets_indices = 1:size(AllData,2);
nlags = 5;
 
[Corrvals, pvals] = corrcoef(AllData,'Rows','pairwise');
  
for t = 0:nlags 
 
drivers_lagged = AllData(1:end-t, drivers_indices); 
targets_current = AllData((t+1):end, targets_indices); 
new_variables = [drivers_lagged targets_current]; 

[corrs pvals] = corrcoef(new_variables,'Rows','pairwise');

Corrvals_lagged(:,:,t+1) = corrs(1:25,1:25);
pvals_lagged(:,:,t+1) = pvals(1:25,1:25);
end


for  j = 1:size(AllData,2)
    
    for i = 1:size(AllData,2) %lagged sources (all variables)
        
        vector=reshape(Corrvals_lagged(j,i,1:nlags+1),[1,nlags+1]);
        vectorp=reshape(pvals_lagged(j,i,1:nlags+1),[1,nlags+1]);
        vector(vectorp > 0.01) = nan;
        numnans = length(find(isnan(vector)));
        
        [maxval maxind]=max(abs(vector));
        if vector(maxind)<0
            maxval = -maxval;
        end
        
    
        MaxCorrmatrix(j,i)=maxval;
        MaxCorrlag(j,i)=maxind-1;
    
    end
    
end



I_dom = ent.I_dom_normbyH;
I_dom(I_dom==0)=nan;
I_dom_lag = ent.I_dom_lag;

for j = 1:size(AllData,2)
    I_dom(j,j)=nan;
    I_dom_lag(j,j)=nan; 
    %omit self links
end

%plot matrix of MI and associated lags
figure(12)
subplot(1,3,1)
imagesc(I_dom(:,1:6),'AlphaData',~isnan(I_dom(:,1:6)))
set(gca,'Ytick',1:25,'Yticklabels',varnames_all)
set(gca,'Xtick',1:25,'Xticklabels',varnames_all)
xtickangle(gca,60)
caxis([0 .6])
colorbar('Location','southoutside')
title('MI for Annual Counts')
ylabel('Sources')
xlabel('Targets')
colormap(copper)

subplot(1,3,2)
imagesc(I_dom_lag(:,1:6),'AlphaData',~isnan(I_dom(:,1:6)))
%set(gca,'Ytick',1:25,'Yticklabels',varnames_all)
set(gca,'Xtick',1:25,'Xticklabels',varnames_all)
xtickangle(gca,60)
caxis([0 5])
colorbar('Location','southoutside')
title('Dominant time lags')
xlabel('Targets')
colormap(copper)


%plot matrix of correlations and associated lags
figure(13)
subplot(1,2,1)
imagesc(MaxCorrmatrix,'AlphaData',~isnan(MaxCorrmatrix))
set(gca,'Ytick',1:25,'Yticklabels',varnames_all)
set(gca,'Xtick',1:25,'Xticklabels',varnames_all)
xtickangle(gca,60)
caxis([-1 1])
colorbar('Location','southoutside')
title('Correlations for Annual Counts')
ylabel('Sources')
xlabel('Targets')

subplot(1,2,2)
imagesc(MaxCorrlag,'AlphaData',~isnan(MaxCorrmatrix))
%set(gca,'Ytick',1:25,'Yticklabels',varnames_all)
set(gca,'Xtick',1:25,'Xticklabels',varnames_all)
xtickangle(gca,60)
caxis([0 5])
colorbar('Location','southoutside')
title('Dominant time lags')
xlabel('Targets')

%comparison of cases where both MI and Corr detected, vs one or the other


for i = 1:size(AllData,2)
    for j=1:size(AllData,2)
        if isnan(MaxCorrmatrix(i,j)) && isnan(I_dom(i,j))
            Match(i,j) = nan;
        elseif ~isnan(MaxCorrmatrix(i,j)) && ~isnan(I_dom(i,j))
            Match(i,j) = 0;            
        elseif ~isnan(MaxCorrmatrix(i,j)) && isnan(I_dom(i,j))
            Match(i,j)=-1;            
        else
            Match(i,j)=1;           
        end
    end
end

num_tot_links = sum(sum(~isnan(Match)))
num_both_links = sum(sum(Match == 0))
num_corr_links = sum(sum(Match == -1))
num_MI_links = sum(sum(Match == 1))

%%
figure(12)
subplot(1,3,3)
imagesc(Match(:,1:6), 'AlphaData',~isnan(Match(:,1:6)))
set(gca,'Xtick',1:25,'Xticklabels',varnames_all)
xtickangle(gca,60)
%set(gca,'Ytick',1:25)
colormap(jet)
colorbar('Location','southoutside')
