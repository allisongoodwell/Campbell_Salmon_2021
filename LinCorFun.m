function [lincor] = LinCorFun(mi,Data,seg)

%Function to compute lagged linear correlations to compare with
%entropy/TIPNet measures (which are computed in EntropyFun)

nvars = mi.nvars;
nTests = mi.nTests;
N = mi.N;
lagvect = mi.lagvect;
bin_scheme=mi.bin_scheme;
Range = mi.Range;
z_opt = mi.ZeroLagOpt;
DomNormOpt = 0;

targets = mi.Targets; %added June 2020, only need these variables as targets of info
ntargs = length(targets);

ndata = size(Data,1);

%find variance (corrolary to entropy)
for i =1:nvars  
    dat = Data(:,i);
    
    dat(isnan(dat))=[]; %omit any nan values to get variance
    VarDat(i) = std(dat)./mean(dat); %coefficient of variance

end

nlags = length(lagvect);

C_lags = zeros(nlags,nvars,ntargs); 
C_inst = zeros(nvars,ntargs);
C_dom = zeros(nvars,ntargs);
C_tau = nan(nvars,ntargs);

%% compute correlation and lagged correlations for all lags and pairs

%first compute lagged corr between all targets and Bonneville, find lag of
%Bonneville max influence

B_maxlag=zeros(1,ntargs);

for Target=1:ntargs
    B = Data(:,1);          %bonneville data
    Y = Data(:,Target);     %target dams
    
    inds = find(isnan(B) | isnan(Y));
    B(inds)=[]; Y(inds)=[];
    
    
    [cvals, lags] = xcorr(B,Y,nlags,'coeff');
    
    
    
    CBon(:,Target) = cvals(1:nlags);
                       

    [val B_maxlag(Target)] = max(CBon(:,Target));
       
end

B_maxlag = nlags - B_maxlag;


for Source = 1:nvars %transmitters
    
    X = Data(:,Source);
 
    for Target = 1:ntargs %receivers
    
       Y = Data(:,Target);
        
       inds = find(isnan(X) | isnan(Y));
       Y(inds)=[]; 
       Xnew = X; Xnew(inds)=[];
        
        cvals = xcorr(Xnew,Y,nlags,'coeff');    
         
        C_lags(:,Source,Target) = cvals(1:nlags);
        
        
        C_inst(Source,Target) = xcorr(Xnew,Y,0,'coeff');       
        [C_dom(Source, Target), C_tau(Source, Target)] = max(C_lags(:,Source,Target));

    end                                                             %receivers
end%transmitters

%find p-values for correlations (maximum correlation values)


for Source = 1:nvars %transmitters
    X = Data(:,Source);
    for Target = 1:ntargs %receivers
        Y = Data(:,Target);


        
    end
end


C_tau = nlags - C_tau;

%% save output products in entropy structure
lincor.C_inst           = C_inst;
lincor.C_lags           = C_lags;
lincor.CBon             = CBon;
lincor.B_maxlag         = B_maxlag;
lincor.variance         = VarDat;
lincor.C_dom            = C_dom;
lincor.C_tau            = C_tau;
end

