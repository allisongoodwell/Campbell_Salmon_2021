%%%%%% TIPNet for CRB Salmon study: Results analysis %%%%%
% Allison Goodwell and Nick Campbell
% June 2020
% This code loads in results files from TIPNet_Salmon.m (10 YEAR CASE -
% DAILY DATA)
% and makes some plots of results

% notes: 
% 1. at beginning, must pick which analysis types to look at for some plots
% 2. using normalized MI and TE, both normalized by min source entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

load('RESULTS_10yr.mat')

ent_spring = ent{1}; %first index is spr/sum/fall case, second is values (1) increments (2)
ent_fall = ent{2}; 

lin_spring = lincor{1}; %lagged linear correlation results (to compare with MI)
lin_fall = lincor{2};

Clags_spring = lin_spring.C_lags;
Clags_fall = lin_fall.C_lags;


Ilags_spring = ent_spring.I_tau_normbyH; %lags, source, target
Ilags_spring(Ilags_spring==0)=nan;

Ilags_fall = ent_fall.I_tau_normbyH; %lags, source, target
Ilags_fall(Ilags_fall==0)=nan;

TElags_spring = ent_spring.TE_Bonneville; %lags, source, target
TElags_spring(TElags_spring==0)=nan;

TElags_fall = ent_fall.TE_Bonneville; %lags, source, target
TElags_fall(TElags_fall==0)=nan;

Bonlag_spring = ent_spring.Bondville_MaxLag;
Bonlag_fall = ent_fall.Bondville_MaxLag;

nsources =8;
ntargets =6;
source_inds = [1 2 3 4 5 6 7 10]; %only look at 1 temperature and one flow variable

colorvect = jet(nsources); 
colorvect(5,:)=[.8 .8 0];

figure(1) %plot sources in each panel
for i = 1:nsources
    si = source_inds(i); %index of source variable
    
    subplot(2,8,i)
    for j = 1:ntargets
        
        if si ~= j
            vects = reshape(Ilags_spring(:,si,j),[31,1]);
            plot(vects,'LineWidth',2,'Color',colorvect(j,:))
        end
        hold on
        vects = reshape(TElags_spring(:,si,j),[31,1]);
        plot(vects,'LineWidth',2,'LineStyle',':','Color',colorvect(j,:))
    end
    title(['spring source ' varnames_all(si)])
    ylim([0 .5])
     xlim([0 30])
    
    subplot(2,8,i+8)
    for j = 1:ntargets
        if si ~= j
            vects = reshape(Ilags_fall(:,si,j),[31,1]);
            plot(vects,'LineWidth',2,'Color',colorvect(j,:))
        end
        hold on
        vects = reshape(TElags_fall(:,si,j),[31,1]);
        plot(vects,'LineWidth',2,'LineStyle',':','Color',colorvect(j,:))
    end
    
    title(['fall source ' varnames_all(si)])
    ylim([0 .5])
    xlim([0 30])
    
end

%% specific figures

%BON as an MI source in fall and spring
figure(2) %plot sources in each panel
ct=1;
for j = 2:ntargets
    
    subplot(2,2,1)
    hold on
    vects = reshape(Ilags_spring(:,1,j),[31,1]);
    hh(j-1)=plot(vects,'LineWidth',2,'Color',colorvect(j,:));
    plot(Bonlag_spring(j),max(vects),'.','Markersize',20,'Color',colorvect(j,:))
    
    subplot(2,2,2)
    hold on
    vects = reshape(Ilags_fall(:,1,j),[31,1]);
    ct=ct+1;
    plot(vects,'LineWidth',2,'Color',colorvect(j,:));
    plot(Bonlag_fall(j),max(vects),'.','Markersize',20,'Color',colorvect(j,:))    
    ct=ct+1;
end

legend(hh,{'MCN','IHR','LWG','PRD','WEL'})
subplot(2,2,1)
title(['MI(BON;salmon) spring'])
ylim([0 .55])
xlim([0 30])
xlabel('time lag (days)')
ylabel('normalized MI (bits/bit)')

subplot(2,2,2)
title(['MI(BON;salmon) fall'])
ylim([0 .55])
xlim([0 30])
xlabel('time lag (days)')
ylabel('normalized MI (bits/bit)')

%correlation version

ct=1;
for j = 2:ntargets

    subplot(2,2,3)
    hold on
    vects = (reshape(Clags_spring(:,1,j),[31,1]));
    hh(j-1)=plot(flipud(vects),'LineWidth',2,'Color',colorvect(j,:));

    subplot(2,2,4)
    hold on
    vects = (reshape(Clags_fall(:,1,j),[31,1]));
    ct=ct+1;
    plot(flipud(vects),'LineWidth',2,'Color',colorvect(j,:));   
    ct=ct+1;
end

legend(hh,{'MCN','IHR','LWG','PRD','WEL'})
subplot(2,2,3)
title(['Corr(BON;salmon) spring'])

xlim([0 30])
xlabel('time lag (days)')
ylabel('corr coeff')

subplot(2,2,4)
title(['Corr(BON;salmon) fall'])

xlim([0 30])
xlabel('time lag (days)')
ylabel('corr coeff')


%%  UPDATE - Q and T (old code below for TE and MI measures)

figure(3)

Vars = categorical(fishsitenames);
Vars = reordercats(Vars,fishsitenames);

subplot(2,3,1)
bar(Vars,[ent_spring.U_T_Bon; ent_spring.U_Bon_T; ent_spring.S_BonT; ent_spring.R_BonT]','stacked')
ylim([0 .8])
title('S_1: T_{lower} and S_2: BON Counts')
legend({'U_{s1}','U_{s2}','S','R'})


subplot(2,3,2)
bar(Vars,[ent_spring.U_Q_Bon; ent_spring.U_Bon_Q; ent_spring.S_BonQ; ent_spring.R_BonQ]','stacked')
ylim([0 .8])
title('S_1: Q_{lower} and S_2: BON Counts')


subplot(2,3,3)
bar(Vars,[ent_spring.U_T_Q; ent_spring.U_Q_T; ent_spring.S_QT; ent_spring.R_QT]','stacked')
ylim([0 .8])
title('S_1: T_{lower} and S_2: Q_{lower}')

subplot(2,3,4)
bar(Vars,[ent_fall.U_T_Bon; ent_fall.U_Bon_T; ent_fall.S_BonT; ent_fall.R_BonT]','stacked')
ylim([0 .8])

subplot(2,3,5)
bar(Vars,[ent_fall.U_Q_Bon; ent_fall.U_Bon_Q; ent_fall.S_BonQ; ent_fall.R_BonQ]','stacked')
ylim([0 .8])

subplot(2,3,6)
bar(Vars,[ent_fall.U_T_Q; ent_fall.U_Q_T; ent_fall.S_QT; ent_fall.R_QT]','stacked')
ylim([0 .8])




%%


%figure 4...influences from neighbor to neighbor
%BON to MCN (1,2)
%MCN to IHR (2,3)
%MCN to PRD (2,5)
%IHR to LWG (3,4)
%PRD to WEL (5,6)
%fall and spring, TE and MI (4 lines per plot)
figure(4)
pairs = [1 2; 2 3; 2 5; 3 4; 5 6];
for i =1:5
    sourcevar = pairs(i,1);
    tarvar = pairs(i,2);
    
    subplot(1,5,i)
        vect = reshape(Ilags_spring(:,sourcevar,tarvar),[31,1]);
    hold on
    plot(vect,'-','Color',colorvect(sourcevar,:),'LineWidth',2)
    vect = reshape(TElags_spring(:,sourcevar,tarvar),[31,1]);
    plot(vect,'x','Color',colorvect(sourcevar,:))
    vect = reshape(Ilags_fall(:,sourcevar,tarvar),[31,1]);
    plot(vect,'--','Color',colorvect(sourcevar,:),'LineWidth',2)
    vect = reshape(TElags_fall(:,sourcevar,tarvar),[31,1]);
    plot(vect,'o','Color',colorvect(sourcevar,:))
    title(sprintf('%s to %s',varnames_all{sourcevar},varnames_all{tarvar}))
    ylim([0 .55]);
    xlim([0 30]);
    if i==5
    legend('MI spring','TE spring','MI fall','TE fall')
    end
end

%% plots of cross-tributary interactions

figure(5)

pairs = [3 5; 3 6; 5 3; 5 4];
for i =1:4
    sourcevar = pairs(i,1);
    tarvar = pairs(i,2);
    
    subplot(2,2,i)
    vect = reshape(Ilags_spring(:,sourcevar,tarvar),[31,1]);
    hold on
    plot(vect,'-','Color',colorvect(sourcevar,:),'LineWidth',2)
    vect = reshape(TElags_spring(:,sourcevar,tarvar),[31,1]);
    plot(vect,'x','Color',colorvect(sourcevar,:))
    vect = reshape(Ilags_fall(:,sourcevar,tarvar),[31,1]);
    plot(vect,'--','Color',colorvect(sourcevar,:),'LineWidth',2)
    vect = reshape(TElags_fall(:,sourcevar,tarvar),[31,1]);
    plot(vect,'o','Color',colorvect(sourcevar,:))
    title(sprintf('%s to %s',varnames_all{sourcevar},varnames_all{tarvar}))
    ylim([0 .55]);
    xlim([0 30]);
    legend('MI spring','TE spring','MI fall','TE fall')
end





%% make some plots of the data...

%15 columns, column 1 is the date
data = xlsread('masterfile10yr_b.xlsx');

date = data(:,1);
salmon_pops = data(:,2:7);
flows = data(:,8:10);
temps = data(:,[14 12 15]); 


%want to convert date to decimal year...
[Y, M, D]=datevec(date-1); %don't know why datevec gives dates 1 day off of excel, -1 to fix
dy = decyear(Y,M,D)+1900; %obtain decimal year (starts 2009)

fishsitenames = {'BON','MCN','IHR','LWG','PRD','WEL'};
flownames = {'LowerQ','SnakeQ','UpperQ'};
tempnames = {'LowerT','SnakeT','UpperT'};

%now want to segment by seasons
%will create 3 datasets: spring, summer/fall, and spring/summer/fall
%breakpoints: DOY 60 (March 1), DOY 180 (late June), DOY 330 (late Nov)
x1=90; x2=220; x3=330;

%need to convert date vector into DOY vector
jd = juliandate([Y+1900,M,D]);

for i = 1:length(date)
    yr = Y(i);
    jd1 = juliandate([yr+1900,1,1]);
    DOY(i) = jd(i)-jd1;
end

%find daily mean, max, min, and 75th and 25th percentiles for each day of
%year
figure(6)
for site =1:6
    
    fish = salmon_pops(:,site);
    
    for d = 1:365
        pop_day = fish(DOY==d);
        avg_fish(d)=median(pop_day);
        max_fish(d)=max(pop_day);
        min_fish(d)=min(pop_day);
        fish75(d) = prctile(pop_day,75);
        fish25(d) = prctile(pop_day,25);
    end
    
    subplot(1,6,site)
    
    hold on
    
    patch([1:365 fliplr(1:365)], [min_fish fliplr(max_fish)],colorvect(site,:))
    plot(1:365,avg_fish,'k','LineWidth',2)
    line([220 220],[0 max(max_fish)],'Linewidth',2','LineStyle','--','Color','k')
    %set(gca,'XTickLabels',[])
    %set(gca,'YTickLabels',[])
    xlim([90 330])
    
end

figure(7)

Q = flows(:,1);
T = temps(:,1);

    for d = 1:365
        q_day = Q(DOY==d);
        t_day = T(DOY==d);
        
        avg_q(d)=median(q_day);
        max_q(d)=max(q_day);
        min_q(d)=min(q_day);
        avg_t(d)=median(t_day);
        max_t(d)=max(t_day);
        min_t(d)=min(t_day);

    end
    
    subplot(1,2,1)   
    patch([1:365 fliplr(1:365)], [min_t fliplr(max_t)],colorvect(7,:))
    hold on
    plot(1:365,avg_t,'k','LineWidth',2)
    line([220 220],[0 max(max_t)],'Linewidth',2','LineStyle','--','Color','k')
    xlim([90 330])
    ylim([40 80])
    title('Temperature')
    
    subplot(1,2,2)   
    patch([1:365 fliplr(1:365)], [min_q fliplr(max_q)],colorvect(8,:))
    hold on
    plot(1:365,avg_q,'k','LineWidth',2)
    line([220 220],[0 max(max_q)],'Linewidth',2','LineStyle','--','Color','k')
    xlim([90 330])
    title('Flow')
    

%% create text files for Circos graphs (uncomment if desired)

%first get matrix of maximum values and associated lags 


for source =1:nsources
    
    si = source_inds(source); %index of source variable
    for target = 1:ntargets
        
        MI_vect_spring = Ilags_spring(:,si,target);
        TE_vect_spring = TElags_spring(:,si,target);
        Corr_vect_spring = Clags_spring(:,si,target);
        
        MI_vect_spring(isnan(MI_vect_spring))=0;
        TE_vect_spring(isnan(TE_vect_spring))=0;
        Corr_vect_spring(isnan(Corr_vect_spring))=0;
        
        
        [MI_matrix_spring(source,target), lagMI_matrix_spring(source,target)] = max(MI_vect_spring);
        [TE_matrix_spring(source,target), lagTE_matrix_spring(source,target)] = max(TE_vect_spring);
        [Corr_matrix_spring(source,target), lagCorr_matrix_spring(source,target)] = max(Corr_vect_spring);
        
        MI_vect_fall = Ilags_fall(:,si,target);
        TE_vect_fall = TElags_fall(:,si,target);
        Corr_vect_fall = Clags_fall(:,si,target);
        
        MI_vect_fall(isnan(MI_vect_fall))=0;
        TE_vect_fall(isnan(TE_vect_fall))=0;
        Corr_vect_fall(isnan(Corr_vect_fall))=0;
        
        
        [MI_matrix_fall(source,target), lagMI_matrix_fall(source,target)] = max(MI_vect_fall);
        [TE_matrix_fall(source,target), lagTE_matrix_fall(source,target)] = max(TE_vect_fall);
        [Corr_matrix_fall(source,target), lagCorr_matrix_fall(source,target)] = max(Corr_vect_fall);
        
    end
end

%set MI self links (=1) to zero
for target=1:ntargets
    MI_matrix_spring(target,target)=0;
    MI_matrix_fall(target,target)=0;
    Corr_matrix_fall(target,target)=0;
    Corr_matrix_spring(target,target)=0;
end


varnames_all = {'BON','MCN','IHR','LWG','PRD','WEL','Q','T'};

%%

foldername = 'CircosFiles';
if exist(foldername)~=7
    mkdir(foldername);
end
addpath(foldername);

colR = round(jet(nsources).*255);
colSource = 255;

for j =1:6
    
    if j==1
        fileID = fopen([foldername,'/MIfall_10yr','.txt'],'w+');
        matrix = MI_matrix_fall;
    elseif j==2
        fileID = fopen([foldername,'/TEfall_10yr','.txt'],'w+');
        matrix = TE_matrix_fall;
    elseif j==3
        fileID = fopen([foldername,'/MI_spring_10yr','.txt'],'w+');
        matrix = MI_matrix_spring;
    elseif j==4
        fileID = fopen([foldername,'/TE_spring_10yr','.txt'],'w+');
        matrix = TE_matrix_spring;
    elseif j==5
        fileID = fopen([foldername,'/Corr_spring_10yr','.txt'],'w+');
        matrix = Corr_matrix_spring;
    elseif j==6
        fileID = fopen([foldername,'/Corr_fall_10yr','.txt'],'w+');
        matrix = Corr_matrix_fall;             
    end
    
    tot_info(j) = sum(sum(matrix));
    
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
    
    for vp = 1:size(matrix,1)
        fprintf(fileID,'%d,%d,%d %s ',colR(vp,1), colR(vp,2), colR(vp,3), varnames_all{vp});

        for v =1:ntargets
            fprintf(fileID,'%d  ',round(matrix(vp,v).*100));
        end

        fprintf(fileID,'\n');
    end


end






