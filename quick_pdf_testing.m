%pdf testing...1D histogram...

clear all
close all
clc

Data = rand(1,10000)';
N=15;
bin_scheme = 'global';
Range = [0; 1];

[pdf, Coords]= compute_pdfGUI(Data,N,bin_scheme,Range);

figure(1)
subplot(2,1,1)
bar(pdf)

subplot(2,1,2)
hist(Data,15)


%%
data = xlsread('masterfile10yr_b.xlsx');

date = data(:,1);
salmon_pops = data(:,2:7);
flows = data(:,8:10);
temps = data(:,[14 12 15]); 

AllData = [log(salmon_pops), log(flows), log(temps)];

AllData(AllData<0)=0;
x1=90; x2=220; x3=330;

%want to convert date to decimal year...
[Y, M, D]=datevec(date-1); %don't know why datevec gives dates 1 day off of excel, -1 to fix
dy = decyear(Y,M,D)+1900; %obtain decimal year (starts 2009)


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

BON_spr = AllData_spr(:,1);
%BON_spr(isnan(BON_spr))=[];

Range = [mi.Range(1,1); mi.Range(2,1)];


[pdf, Coords]= compute_pdfGUI(BON_spr,N,bin_scheme,Range);

figure(2)
subplot(2,1,1)
bar(pdf)

subplot(2,1,2)
hist(BON_spr,15)
