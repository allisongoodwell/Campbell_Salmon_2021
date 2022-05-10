function [entropy] = EntropyFun(mi,Data,seg)
%Function to compute various entropy, transfer entropy, information decomposition measures
%Allison Goodwell, June 2014

%Update 7/1/20: altered code for Nick Salmon cases,
%and adding TE_Bonneville variable

%Update 9/21/20: alter code for Salmon, change TE_Bonneville to condition
%on dominant lag for Bonneville to given source....

nvars = mi.nvars;
nTests = mi.nTests;
N = mi.N;
lagvect = mi.lagvect;
bin_scheme=mi.bin_scheme;
Range = mi.Range;
z_opt = mi.ZeroLagOpt;
DomNormOpt = 0;
z_effect = mi.DataPrep.Z_effect;

targets = mi.Targets; %added June 2020, only need these variables as targets of info
ntargs = length(targets);

ndata = size(Data,1);

if z_opt==1 %include zero in lagvect
lagvect =[0 lagvect];
end

nlags = length(lagvect);
TE = zeros(nlags,nvars,ntargs);
TE_normbyItot = TE;
I_lags = zeros(nlags,nvars,ntargs);
I_lags_nonsig = zeros(nlags,nvars,ntargs);
I_normbyH = zeros(nlags,nvars,ntargs);
I_inst = zeros(nvars,ntargs);
I_inst_normbyH = zeros(nvars,ntargs);
I_inst_sig = zeros(nvars,ntargs);
I_pvalue = zeros(nvars,ntargs);
H_x_1 = zeros(nvars,1);
H_x_2 = zeros(nvars,1);

TE_T_normbyItot = zeros(nvars,ntargs);


TI = zeros(nvars,ntargs);
TE_T = TI;
I_dom = zeros(nvars,ntargs);
I_dom_normbyH = zeros(nvars,ntargs);
I_tau = nan(nvars,ntargs);

U_Svar_Bon = zeros(nvars,ntargs);
U_Bon_Svar = zeros(nvars,ntargs);
R_BonSvar = zeros(nvars,ntargs);
S_BonSvar = zeros(nvars,ntargs);

%AEG: adding a variable for TE --> conditions on the current value of the first
%node, instead of past value of target node
%also adding lag that corresponds to max value
TE_Bonneville = zeros(nlags,nvars,ntargs);


%% compute 1D entropy values for each node
for Source = 1:nvars
    X = Data(:,Source);
    pdf = compute_pdfGUI(X,N,bin_scheme, Range(:,Source));
    info = compute_info_measures(pdf);
    H_x_1(Source) = info.Hx;
    H_pdf_vals(:,Source) = pdf; %save 1D pdfs for each variable (N bins)
end


%% compute MI and lagged MI for all lags and pairs

timer=timebar(sprintf('MI for segment %d',seg),'Network Timer');
ct=0;
timebar(timer,ct/nvars^2);

%update 9/2020: first compute lagged MI between all targets and Bonneville, find lag of
%Bonneville max influence

B_maxlag=zeros(1,ntargs);
for Target=1:ntargs
    B = Data(:,1); %bonneville data
    Y = Data(:,Target); %target dams
    
    for t=1:nlags %may or may not include zero lag           

            lagt = lagvect(t);

            nTuples = ndata-lagt-1;
            tar_start=1+lagt;
            
            Svar=B(1:nTuples);                              %Leading Node S (lag tau earlier than present)
            Tar=Y(tar_start:tar_start+nTuples-1);           %Led Node Tar (one timestep in future)
            
            Tuple = [Svar Tar];
            Tuple(all(isnan(Tuple),2),:) = [];
        
            pdf = compute_pdfGUI(Tuple,N,bin_scheme, [Range(:,1)  Range(:,Target)]);
            info = compute_info_measures(pdf);
            I_lags_BON(t)=info.I;
                                  
    end
    
    [~, B_maxlag(Target)] = max(I_lags_BON); %identify lag at which BON most informs all salmon counts
       
end

B = Data(:,1); %Bonneville salmon data is first index 
%mutual information between sources and targets
for Source = 1:nvars %transmitters
        
    X = Data(:,Source); 
        
    for Target = 1:ntargs %receivers

        BonLag = B_maxlag(Target);    
        Y = Data(:,Target);
        
        U_Svar_Bon_vect = zeros(nlags); %unique from source for lag
        U_Bon_Svar_vect = zeros(nlags);  %unique from BON for lag
        R_BonSvar_vect = zeros(nlags); 
        S_BonSvar_vect = zeros(nlags); 
        Itot_vect = zeros(nlags); 

        for t=1:nlags %may or may not include zero lag           

            lagt = lagvect(t);
            
            %first compute MI
            nTuples = ndata-lagt-1;
            tarlag_start=1;
            tar_start=1+lagt;
                      
            Svar=X(1:nTuples);                              %Leading Node S (lag tau earlier than present)
            Tar=Y(tar_start:tar_start+nTuples-1);           %Led Node Tar (one timestep in future)
     
            Tuple = [Svar Tar];
            Tuple(all(isnan(Tuple),2),:) = [];
     
            pdf = compute_pdfGUI(Tuple,N,bin_scheme, [Range(:,Source)  Range(:,Target)]);
            info = compute_info_measures(pdf);
            I_lags(t,Source,Target)=info.I;
            Htemp = min(info.Hx1,info.Hx2);
            
            I_lags_nonsig(t,Source,Target) = info.I; %retain all values for later S-U-R...
                    
            %shuffle to compute significance
            if nTests>0
                if t==1 %do shuffled significance testing for zero lag only
                    I_shuff=zeros(1,nTests);
                    for test=1:nTests
                        Tar_shuff = randsample(Tar,length(Tar));
                        
                        Tuple = [Svar Tar_shuff];
                        Tuple(all(isnan(Tuple),2),:) = [];
                        
                        pdfshuff = compute_pdfGUI(Tuple,N,bin_scheme, ...
                            [Range(:,Source)  Range(:,Target)]);
                        infoshuff = compute_info_measures(pdfshuff);
                        I_shuff(test) = infoshuff.I;
                    end
                    I_shuff_sig = mean(I_shuff)+3*std(I_shuff);
                    [h,p,ci,stats] = ttest((info.I-I_shuff)./std(I_shuff));
                end

                I_pvalue(Source,Target)=p;
                I_inst_sig(Source,Target)=I_shuff_sig;
 
                if info.I-I_shuff_sig< 0
                I_lags(t,Source,Target) = 0; %set to zero if non-stat-sig link
                end    
            end
            
            I_normbyH(t,Source,Target)=I_lags(t,Source,Target)./Htemp;
             
            %next compute TE, conditioned on BON
            nTuples = ndata-max(lagt,BonLag)-1;
                        
            if lagt>BonLag
                tarlag_start=1;
                Bon_start=1+lagt-BonLag;
                tar_start=1+lagt;
            else
                Bon_start =1;
                tarlag_start = 1+BonLag-lagt;
                tar_start=1+BonLag;
            end
            

            Bonn=B(Bon_start:Bon_start+nTuples-1);          %AEG added: Led Node B (Bonneville)
            Svar=X(1:nTuples);                              %Leading Node S (lag tau earlier than present)
            Tar=Y(tar_start:tar_start+nTuples-1);           %Led Node Tar (one timestep in future)
            
            
            % for significant lagged I, also compute Transfer Entropy TE
            if I_lags(t,Source,Target)>0 && t>=BonLag
                               
                %TE_Bonneville Calculation here
                Tuple = [Svar Bonn Tar];
                Tuple(all(isnan(Tuple),2),:) = [];
                
                pdf = compute_pdfGUI(Tuple,N,bin_scheme,...
                    [Range(:,Source) Range(:,1) Range(:,Target)]);
                info = compute_info_measures(pdf);
                
                Htemp = min(info.Hx1,info.Hx2);
                
                TE_Bonneville(t,Source,Target)=info.T./Htemp;
                
                %information decomposition               
                U_Svar_Bon_vect(t) = info.U1/info.Hy; %unique from source for lag
                U_Bon_Svar_vect(t) = info.U2/info.Hy; %unique from BON for lag
                R_BonSvar_vect(t) = info.R/info.Hy;
                S_BonSvar_vect(t)=info.S/info.Hy;
                Itot_vect(t) = info.Itot/info.Hy;
               
            end
                   
        end                                                         %lags
        
        I_inst(Source,Target) = I_lags(1,Source,Target);
        I_inst_normbyH(Source,Target) = I_normbyH(1,Source,Target);
      
        %find lag of max Itot for each source,target pair
        lagind = find(Itot_vect == max(Itot_vect));
        
        if size(lagind)>0
            lagind = lagind(1);
        
            MaxLag(Source,Target) = lagind;
        
            U_Svar_Bon(Source,Target) = U_Svar_Bon_vect(lagind);
            U_Bon_Svar(Source,Target) = U_Bon_Svar_vect(lagind);
            R_BonSvar(Source,Target) = R_BonSvar_vect(lagind);
            S_BonSvar(Source,Target) = S_BonSvar_vect(lagind);
        else
            MaxLag(Source,Target) = nan;
        
            U_Svar_Bon(Source,Target) = 0;
            U_Bon_Svar(Source,Target) = 0;
            R_BonSvar(Source,Target) = 0;
            S_BonSvar(Source,Target) = 0;
        end
        
        
        ct=ct+1;
        timebar(timer,ct/nvars^2)
    end                                                             %receivers
end                                                                 %transmitters

%also do SUR for I(Temp, Q; Target) for particular combo of flow rate and
%temperature
Q_ind = 7;
T_ind = 10; 

for Target = 1:ntargs
    
    Tdata = Data(:,T_ind);
    Qdata = Data(:,Q_ind);
    Y = Data(:,Target);
        
    laggedI_q = I_lags_nonsig(:,Q_ind,Target);
    laggedI_t = I_lags_nonsig(:,T_ind,Target);
        
    if max(laggedI_q)>0 && max(laggedI_t)>0
        lagQ = find(laggedI_q==max(laggedI_q));
        lagT = find(laggedI_t==max(laggedI_t));
        nTuples = ndata-max(lagT,lagQ)-1;
        
        if lagT>lagQ       
            T_start=1;
            Q_start=1+lagT-lagQ;
            tar_start=1+lagT;
            
        else
            Q_start =1;
            T_start = 1+lagQ-lagT;
            tar_start=1+lagQ;
        end
               
        Tvar=Tdata(T_start:T_start+nTuples-1);
        Qvar=Qdata(Q_start:Q_start+nTuples-1);
        Tar=Y(tar_start:tar_start+nTuples-1);
        
        Tuple = [Tvar Qvar Tar];
        Tuple(all(isnan(Tuple),2),:) = [];
        
        pdf = compute_pdfGUI(Tuple,N,bin_scheme,...
            [Range(:,T_ind) Range(:,Q_ind) Range(:,Target)]);
        info = compute_info_measures(pdf);
          
        %add SUR to get U_svar, U_bon, R and S
        U_T_Q(Target) = info.U1/info.Hy;
        U_Q_T(Target) = info.U2/info.Hy;
        R_TQ(Target) = info.R/info.Hy;
        S_TQ(Target)=info.S/info.Hy;
             
    end
    
end

close(timer)

%% save output products in entropy structure
entropy.H_x_1 =             H_x_1;              %entropy based on 1D pdf
entropy.H_x_2 =             H_x_2;              %entropy based on 2D pdf
entropy.I_tau =             I_lags;             %mutual info I(lag,source,target)
entropy.I_tau_normbyH =     I_normbyH;          %normalized mutual info I(lag,source,target)
entropy.I_dom =             I_dom;              %maximum lagged mutual information (bits)
entropy.I_dom_normbyH =     I_dom_normbyH;      %maximum lagged mutual information (normalized by H)
entropy.I_dom_lag=          I_tau;              %lag associated with I_dom (timesteps)

entropy.I_inst =            I_inst;             %lag zero mutual information
entropy.I_inst_normbyH=     I_inst_normbyH;     %normalized lag zero mutual information
entropy.I_inst_sigthresh =  I_inst_sig;         %statistical significance threshold for mutual information

entropy.TE_Bonneville =     TE_Bonneville;      %AEG added for Nick salmon thesis
entropy.Bondville_MaxLag =  B_maxlag; 

entropy.H_pdf_vals    =     H_pdf_vals;

entropy.USvarBon =          U_Svar_Bon;

entropy.U_T_Bon =           U_Svar_Bon(T_ind,:);
entropy.U_Bon_T =           U_Bon_Svar(T_ind,:);
entropy.S_BonT =            S_BonSvar(T_ind,:);
entropy.R_BonT =            R_BonSvar(T_ind,:);

entropy.U_Q_Bon =           U_Svar_Bon(Q_ind,:);
entropy.U_Bon_Q =           U_Bon_Svar(Q_ind,:);
entropy.S_BonQ =            S_BonSvar(Q_ind,:);
entropy.R_BonQ =            R_BonSvar(Q_ind,:);

entropy.U_Q_T = U_Q_T;
entropy.U_T_Q = U_T_Q;
entropy.S_QT = S_TQ;
entropy.R_QT = R_TQ;

entropy.I_lags_nonsig = I_lags_nonsig;
entropy.I_pvalue = I_pvalue;

end

