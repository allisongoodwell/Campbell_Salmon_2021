function [entropy] = EntropyFun(mi,Data,seg)
%Function to compute various entropy, transfer entropy, information decomposition measures
%Allison Goodwell, June 2014
%August 2014: also calculate lagged mutual information
%Jan 2015: updated method of KDE
%Jan 2015: alter method of transfer entropy: compute T based on
%conditioning lag of min[H(Yt)|H(Yt-lag))]
%Feb 2015: bin_scheme: 'local' or 'global', Range: min and max data values
%Feb 2015: changed kd combined function to normalize variables
%June 2015: changed KDE and info measures code (to also compute
%redundant/synergistic information
%July 2015: save certain computed pdfs according to [trans, rec, lag]
%vector
%July 2015: consider multiple sources as contributors to target
%maxS = number of sources to any single target varible to compute unique
%information
%August 2015: modify code to reduce number of computations to redundancy:
%first eliminate redundant links from each source node (with multiple
%signifcant lags - many lags likely to be completely redundant)
%September 2015: added cumulative synergy, changed outputs and inputs to
%structures: mi has all parameters, data is matrix of data, entropy
%is structure with entropy results and pdfs according to what should be
%saved
%September 8, 2015: implement alternate version of redundancy (also look
%for correlation between sources)
%September 30, 2015: updated compute_info_measures for output as structure
%January 18, 2016: updates for GUI - input h (smoothing parameter) to KDE
%computation, option for method in mi - fixed binning or KDE
%February 2016: update options: mi.netopt (1=full network, 2 = H and I
%only, 3 = H only (only entropy of nodes)
%February 2016: alter redundancy: 2 redundancy matrices: R_tau is
%redundancy between a single source to a target node at different time
%lags, R_T is redundancy between dominant link from each source to a target

%Update: altered version to consider zero-lag I as possible dominant link
%for T/I, other measures, depending on mi.ZeroLagOpt

%Update 6/10/16: mi.DomNormOpt = 0 for classifying dominant link as
%non-normalized (default), = 1 for normalized (previous method)
%changed several variable names

%Update 7/1/20: altered code for Nick Salmon cases, omitting SUR code
%and adding TE_Bonneville variable

%Update 9/21/20: alter code for Salmon, change TE_Bonneville to condition
%on dominant lag for Bonneville to given source....

nvars = mi.nvars;
nTests = mi.nTests;
N = mi.N;
lagvect = mi.lagvect;
bin_scheme=mi.bin_scheme;
Range = mi.Range;
method = mi.method;
z_opt = mi.ZeroLagOpt;
DomNormOpt = 0;
z_effect = mi.DataPrep.Z_effect;

targets = mi.Targets; %added June 2020, only need these variables as targets of info
ntargs = length(targets);

ndata = size(Data,1);

if z_opt==1 %include zero in lagvect
lagvect =[0 lagvect];
end

for i =1:nvars
  
    dat = Data(:,i);
    Ra = range(dat);
    varterm = var(dat);
    varterm(isnan(varterm))=0;
    
h1D(i)= 1.06 * ndata^(-1/5) .* varterm./Ra; 
h2D(i)= 1.77 * ndata^(-1/6) .* varterm./Ra; 
h3D(i) = 2.78 * ndata^(-1/7) .* varterm./Ra; 
end

nlags = length(lagvect);
lim = eps*10; % machine limit of matlab 
N_sources = zeros(nvars,1);
N_sourcenodes = zeros(nvars,1);
TE = zeros(nlags,nvars,ntargs);
TE_normbyItot = TE;
I_lags = zeros(nlags,nvars,ntargs);
I_normbyH = zeros(nlags,nvars,ntargs);
I_inst = zeros(nvars,ntargs);
I_inst_normbyH = zeros(nvars,ntargs);
I_inst_sig = zeros(nvars,ntargs);
H_x_1 = zeros(nvars,1);
H_x_2 = zeros(nvars,1);

TE_T_normbyItot = zeros(nvars,ntargs);


TI = zeros(nvars,ntargs);
TE_T = TI;
I_dom = zeros(nvars,ntargs);
I_dom_normbyH = zeros(nvars,ntargs);
I_tau = nan(nvars,ntargs);

%AEG: adding a variable for TE --> conditions on the current value of the first
%node, instead of past value of target node
%also adding lag that corresponds to max value
TE_Bonneville = zeros(nlags,nvars,ntargs);


%% compute 1D entropy values for each node
for Source = 1:nvars
    
    X = Data(:,Source);
    
    pdf = compute_pdfGUI(X,N,bin_scheme, Range(:,Source),method,z_effect(Source+1));
    info = compute_info_measures(pdf);
    H_x_1(Source) = info.Hx;
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
            
            Svar=B(1:nTuples);        %Leading Node S (lag tau earlier than present)
            Tar=Y(tar_start:tar_start+nTuples-1);        %Led Node Tar (one timestep in future)
            
            
            Tuple = [Svar Tar];
            Tuple(all(isnan(Tuple),2),:) = [];
        
            pdf = compute_pdfGUI(Tuple,N,bin_scheme, [Range(:,Source)  Range(:,Target)],method,[z_effect(Source+1) z_effect(Target+1)]);
            info = compute_info_measures(pdf);
            I_lags_prelim(t)=info.I;
            
                        
    end
    
    [val B_maxlag(Target)] = max(I_lags_prelim);
       
end






for Source = 1:nvars %transmitters
    
    %AEG: adding first index (assume Bonneville) as B
    B = Data(:,1); 
    
    
    X = Data(:,Source);
    
    
    for Target = 1:ntargs %receivers

        BonLag = B_maxlag(Target);
    
        Y = Data(:,Target);

        for t=1:nlags %may or may not include zero lag           

            lagt = lagvect(t);

            
            %first compute MI
            nTuples = ndata-lagt-1;
            tarlag_start=1;
            tar_start=1+lagt;
            
            
            Svar=X(1:nTuples);        %Leading Node S (lag tau earlier than present)
            Tar=Y(tar_start:tar_start+nTuples-1);        %Led Node Tar (one timestep in future)
     
            Tuple = [Svar Tar];
            Tuple(all(isnan(Tuple),2),:) = [];
     
            pdf = compute_pdfGUI(Tuple,N,bin_scheme, [Range(:,Source)  Range(:,Target)],method,[z_effect(Source+1) z_effect(Target+1)]);
            info = compute_info_measures(pdf);
            I_lags(t,Source,Target)=info.I;
            Htemp = min(info.Hx1,info.Hx2);
                    
            %shuffle to compute significance
            if nTests>0
                if t==1 %do shuffled significance testing for zero lag only
                    I_shuff=zeros(1,nTests);
                    for test=1:nTests
                        Tar_shuff = randsample(Tar,length(Tar));
                        
                        Tuple = [Svar Tar_shuff];
                        Tuple(all(isnan(Tuple),2),:) = [];
                        
                        pdfshuff = compute_pdfGUI(Tuple,N,bin_scheme, ...
                            [Range(:,Source)  Range(:,Target)],method,[z_effect(Source+1) z_effect(Target+1)]);
                        infoshuff = compute_info_measures(pdfshuff);
                        I_shuff(test) = infoshuff.I;
                    end
                    I_shuff_sig = mean(I_shuff)+4*std(I_shuff);
                    [h,p,ci,stats] = ttest((info.I-I_shuff)./std(I_shuff));
                end

                I_pvalue(Source,Target)=p;
                I_inst_sig(Source,Target)=I_shuff_sig;
 
                if info.I-I_shuff_sig< 0
                I_lags(t,Source,Target) = 0;
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
            

            Bonn=B(Bon_start:Bon_start+nTuples-1);        %AEG added: Led Node B (Bonneville)
            Svar=X(1:nTuples);        %Leading Node S (lag tau earlier than present)
            Tar=Y(tar_start:tar_start+nTuples-1);        %Led Node Tar (one timestep in future)
            
            
            % for significant lagged I, also compute Transfer Entropy TE
            if I_lags(t,Source,Target)>0 && t>=BonLag
                               
                %TE_Bonneville Calculation here
                Tuple = [Svar Bonn Tar];
                Tuple(all(isnan(Tuple),2),:) = [];
                
                pdf = compute_pdfGUI(Tuple,N,bin_scheme,...
                    [Range(:,Source) Range(:,1) Range(:,Target)],method);
                info = compute_info_measures(pdf);
                
                Htemp = min(info.Hx1,info.Hx2);
                
                TE_Bonneville(t,Source,Target)=info.T./Htemp;
               
            end
                   
        end                                                         %lags
        
        I_inst(Source,Target) = I_lags(1,Source,Target);
        I_inst_normbyH(Source,Target) = I_normbyH(1,Source,Target);
      
        ct=ct+1;
        timebar(timer,ct/nvars^2)
    end                                                             %receivers
end                                                                 %transmitters


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
entropy.I_inst_sigthresh =  I_inst_sig;   %statistical significance threshold for mutual information

entropy.TE_Bonneville =     TE_Bonneville;  %AEG added for Nick salmon thesis
entropy.Bondville_MaxLag =  B_maxlag; 

end

