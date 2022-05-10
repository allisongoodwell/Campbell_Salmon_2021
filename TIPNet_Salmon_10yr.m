%%%%%% TIPNet for CRB Salmon study %%%%%
% Allison Goodwell and Nick Campbell
% June 2020
% This code (1) loads daily data for salmon pop, flow, water temps at CRB sites
% (2) separates data into spring and summer/fall seasons
% (3) computes IT measures: information flows to salmon populations
% (4) computes linear lagged correlations for comparison with IT measures

% notes: 
% 
% two sub-cases currently: one with raw values of data, one with
% increments
% all variables considered as sources, only salmon counts considered as targets
% results from 6 TIPNet runs saved in results mat file, to be loaded into
% analysis program (AnalyzeResults code) to make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

%15 columns, column 1 is the date
data = xlsread('masterfile10yr_b.xlsx');

date = data(:,1);
%want to convert date to decimal year...
[Y, M, D]=datevec(date-1); %don't know why datevec gives dates 1 day off of excel, -1 to fix
dy = decyear(Y,M,D)+1900; %obtain decimal year (starts 2009)

salmon_pops = data(:,2:7);
flows = data(:,8:10);
temps = data(:,[14 12 15]); 



AllData = [log(salmon_pops), flows, temps];

AllData(AllData<0)=0;

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

mi.Range = [min(AllData);max(AllData)];

AllData_spr = AllData;
AllData_spr(DOY<x1 | DOY>x2,:)=nan;

AllData_fall = AllData;
AllData_fall(DOY<x2 | DOY>x3,:)=nan;


%% run TIPnet for different seasons, and linear lagged correlations

for i =1:2 %1 = spr, 2 = sum+fall
    
    if i ==1
      Data = AllData_spr; 
    elseif i==2
      Data = AllData_fall; 
    end
      
    mi.Targets = 1:6; %only targets of interest are the salmon pops...will reduce TIPNet calc time
    
    mi.nvars=12;
    mi.nTests=100;
    mi.N = 11;
    mi.lagvect = 0:30;
    mi.bin_scheme = 'global';   %changed from earlier version, pdf binning over max and min for both seasons
    mi.ZeroLagOpt = 0;
    mi.DataPrep.Z_effect = zeros(1,mi.nvars+1);
    mi.NoSelfOpt=0;
    
    ent{i} = EntropyFun(mi,Data,i);
    lincor{i} = LinCorFun(mi,Data,i);

end

save('RESULTS_10yr.mat','ent','varnames_all','fishsitenames','flownames','tempnames','lincor')






