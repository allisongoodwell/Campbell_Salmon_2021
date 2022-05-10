%%%%%% TIPNet for CRB Salmon study, 35-year cases %%%%%
% Allison Goodwell and Nick Campbell
% June 2020
% This code (1) loads annual  data for salmon pop and other variables in
% CRB
% (3) computes IT measures: information flows to salmon populations and
% between variables

% notes: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

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

%AEG edit 10/5/20: keep BON, but for others - fraction of BON
for i =2:6
   
    vect=salmon_pops(:,i)./salmon_pops(:,1); %

    %omit outliers
    minval = mean(vect)-1.5*iqr(vect);
    maxval = mean(vect)+1.5*iqr(vect);
    vect(vect>maxval)=maxval;
    vect(vect<minval)=minval;
    salmon_pops(:,i)=vect;
    
end


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

varnames_all = [fishsitenames, flownames, tempnames,... 
    PDO_name, PPT_name, Airtemp_name,DegDay_name, SWE_name];

%preliminary data checking
figure(1)
subplot(3,1,1)
plot(year, salmon_pops)
legend(fishsitenames)
subplot(3,1,2)
plot(year,flows)
legend(flownames)
subplot(3,1,3)
plot(year, watertempdays)
legend(tempnames)


%Running TIPNet...
mi.Targets = 1:25; %1:6 to reduce tipnet calcs...

mi.nvars=25;
mi.nTests=500;
mi.N = 3; %can change this...should be between 3-5 bins for pdfs..
mi.lagvect = 0:5; %max lag of 7 years
mi.bin_scheme = 'local';   
mi.method = 'fixed';
mi.ZeroLagOpt = 0;
mi.NoSelfOpt=0;

mi.Range = [min(AllData); max(AllData)];
ent = EntropyFun35yr(mi,AllData,1);
lincor = LinCorFun(mi,AllData,1);

save('RESULTS_35yr.mat','ent','varnames_all','fishsitenames','lincor');






