%%%%%% TIPNet for CRB Salmon study %%%%%
% Allison Goodwell and Nick Campbell
% June 2020
% This code (1) loads daily data for salmon pop, flow, water temps at CRB sites
% (2) separates data into spring and summer/fall seasons
% (3) computes IT measures: information flows to salmon populations

% notes: 
% 1. may want to change break dates between seasons...
% 2. two sub-cases currently: one with raw values of data, one with
% increments
% 3. all variables considered as sources, only salmon considered as targets
% 4. results from 6 TIPNet runs saved in results mat file, to be loaded into
% analysis program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

%15 columns, column 1 is the date
data = xlsread('masterfile10yr_b.xlsx');

date = data(:,1);
salmon_pops = data(:,2:7);
flows = data(:,8:10);
temps = data(:,[14 12 15]); 

AllData = [salmon_pops, flows, temps];

%want to convert date to decimal year...
[Y, M, D]=datevec(date-1); %don't know why datevec gives dates 1 day off of excel, -1 to fix
dy = decyear(Y,M,D)+1900; %obtain decimal year (starts 2009)

fishsitenames = {'BON','MCN','IHR','LWG','PRD','WEL'};
flownames = {'LowerQ','SnakeQ','UpperQ'};
tempnames = {'LowerT','SnakeT','UpperT'};

varnames_all = [fishsitenames, flownames, tempnames];

%preliminary data checking (looks good)
figure(1)
subplot(3,1,1)
plot(dy,AllData(:,1:6))
legend(fishsitenames)
subplot(3,1,2)
plot(dy,AllData(:,7:9))
legend(flownames)
subplot(3,1,3)
plot(dy,AllData(:,10:12))
legend(tempnames)


%now want to segment by seasons
%will create 3 datasets: spring, summer/fall, and spring/summer/fall
%breakpoints: DOY 60 (March 1), DOY 180 (late June), DOY 330 (late Nov)
%x1 = 60; x2 = 180; x3 =330;
x1=90; x2=220; x3=330;

%need to convert date vector into DOY vector
jd = juliandate([Y+1900,M,D]);

for i = 1:length(date)
    yr = Y(i);
    jd1 = juliandate([yr+1900,1,1]);
    DOY(i) = jd(i)-jd1;
end

AllData_spr = AllData;
AllData_spr(DOY<x1 | DOY>x2,:)=nan;

AllData_sumfall = AllData;
AllData_sumfall(DOY<x2 | DOY>x3,:)=nan;

AllData_sprsumfall = AllData;
AllData_sprsumfall(DOY<x1 | DOY>x3,:)=nan;

%Secondary case changing Q and T to be difference from Lower C
Qvals = AllData(:,7:9);
Tvals = AllData(:,10:12);
Qvals(:,2) = abs(Qvals(:,2)-Qvals(:,1));
Qvals(:,3)= abs(Qvals(:,3)-Qvals(:,1));
Tvals(:,2) = abs(Tvals(:,2)-Tvals(:,1));
Tvals(:,3)= abs(Tvals(:,3)-Tvals(:,1));

AllData2 = [salmon_pops, Qvals, Tvals];

AllData_spr2 = AllData2;
AllData_spr2(DOY<x1 | DOY>x2,:)=nan;

AllData_sumfall2 = AllData2;
AllData_sumfall2(DOY<x2 | DOY>x3,:)=nan;

AllData_sprsumfall2 = AllData2;
AllData_sprsumfall2(DOY<x1 | DOY>x3,:)=nan;

%%
%now have three datasets for different seasons/splits, ready for TIPNet...

for i =1:3 %spr, sum+fall, spr+sum+fall, in each case: values vs. increments
    
    if i ==1
      Data1 = AllData_spr; 
      Data2 = AllData_spr2;
    elseif i==2
      Data1 = AllData_sumfall; 
      Data2 = AllData_sumfall2;
    elseif i==3
      Data1 = AllData_sprsumfall; 
      Data2 = AllData_sprsumfall2;
    end
    
    
    mi.Targets = 1:6; %only targets of interest are the salmon pops...can reduce TIPNet calcs
    
    mi.nvars=12;
    mi.nTests=100;
    mi.N = 15;
    mi.lagvect = 0:30;
    mi.bin_scheme = 'local';   
    mi.method = 'fixed';
    mi.ZeroLagOpt = 0;
    mi.DataPrep.Z_effect = zeros(1,mi.nvars+1);
    mi.NoSelfOpt=0;
    
    mi.Range = [min(Data1); max(Data1)];
    ent{i,1} = EntropyFun(mi,Data1,1);
    
    mi.Range = [min(Data2); max(Data2)];
    ent{i,2} = EntropyFun(mi,Data2,1);
    
end

save('RESULTS_10yr.mat','ent','varnames_all','fishsitenames','flownames','tempnames')






