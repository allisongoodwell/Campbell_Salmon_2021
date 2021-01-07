function [entropy] = EntropyFun35yr(mi,Data,seg)
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

for Source = 1:nvars %transmitters
    
    %AEG: adding first index (assume Bonneville) as B
    B = Data(:,1); 
    
    X = Data(:,Source);
    
    if flag(Source)==1
        ct=ct+nvars;
        continue
    end
    
    for Target = 1:ntargs %receivers
        %fprintf('trans = %d, rec = %d\n',Source,Target)
        
        Y = Data(:,Target);
        if flag(Target)==1
            ct=ct+1;
            continue
        end
   
        for t=1:nlags %may or may not include zero lag           

            lagt = lagvect(t);

            nTuples = ndata-lagt-1;
            tar_start=1+lagt;
            tarlag_start = max(1,lagt);
            
            Svar=X(1:nTuples);        %Leading Node S (lag tau earlier than present)
            Tar=Y(tar_start:tar_start+nTuples-1);        %Led Node Tar (one timestep in future)
            TarLag=Y(tarlag_start:tarlag_start+nTuples-1);        %current node TarLag
            Bonn=B(tar_start:tar_start+nTuples-1);        %AEG added: Led Node B (Bonneville)
            
            Tuple = [Svar Tar];
            Tuple(all(isnan(Tuple),2),:) = [];

           
            %first compute mutual information for all lags
            %fprintf('computing lagged I Source %d Target %d lag %d\n',Source,Target,lagt)
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
            
            % for significant lagged I, also compute Transfer Entropy TE
            if I_lags(t,Source,Target)>0
                
                Tuple = [Svar TarLag Tar];
                Tuple(all(isnan(Tuple),2),:) = [];
                
                
                pdf = compute_pdfGUI(Tuple,N,bin_scheme,...
                    [Range(:,Source) Range(:,Target) Range(:,Target)],method,[z_effect(Source+1) z_effect(Target+1)]);
                info = compute_info_measures(pdf);
                TE(t,Source,Target)=info.T;
                if info.Itot>0
                TE_normbyItot(t,Source,Target) = info.T/info.Itot; %normalized by total information shared
                end
                
                %AEG: adding TE_Bonneville Calculation here
                Tuple = [Svar Bonn Tar];
                Tuple(all(isnan(Tuple),2),:) = [];
                
                pdf = compute_pdfGUI(Tuple,N,bin_scheme,...
                    [Range(:,Source) Range(:,1) Range(:,Target)],method);
                info = compute_info_measures(pdf);
                TE_Bonneville(t,Source,Target)=info.T;
               
            end
                   
        end                                                         %lags
        
        if z_opt ==1 %set I_inst to first index of I_lags
            I_inst(Source,Target) = I_lags(1,Source,Target);
            I_inst_normbyH(Source,Target) = I_normbyH(1,Source,Target);
        else   %otherwise compute instantaneous mutual info separately
            
            Tuple = [X Y];
            Tuple(all(isnan(Tuple),2),:) = [];
            
            pdf = compute_pdfGUI(Tuple,N,bin_scheme, [Range(:,Source)  Range(:,Target)],method, [z_effect(Source+1) z_effect(Target+1)]);
            info = compute_info_measures(pdf);
            I_inst(Source,Target)=info.I;
            
            if nTests>0
                I_shuff=zeros(1,nTests);
                for test=1:nTests
                    Tar_shuff = randsample(Y,length(Y));
                    
                    Tuple = [X Tar_shuff];
                    Tuple(all(isnan(Tuple),2),:) = [];
                    
                    pdfshuff = compute_pdfGUI(Tuple,N,bin_scheme, ...
                        [Range(:,Source)  Range(:,Target)],method, [z_effect(Source+1) z_effect(Target+1)]);
                    infoshuff = compute_info_measures(pdfshuff);
                    I_shuff(test) = infoshuff.I;
                end
                I_shuff_sig = mean(I_shuff)+4*std(I_shuff);
                
                I_inst_sig(Source,Target)=I_shuff_sig;
               
                if I_shuff_sig - info.I > 0
                
                I_inst(Source,Target) =0; 
                
                else
                I_inst_normbyH(Source,Target) = I_inst(Source,Target)./min(info.Hx1,info.Hx2);
                end
                
                
            end
            
        end
        
        ct=ct+1;
        timebar(timer,ct/nvars^2)
    end                                                             %receivers
end                                                                 %transmitters



%re-define H_x as I(X;X) (from the 2D pdf instead of 1D)
for n =1:ntargs
H_x_2(n) = I_inst(n,n);
end

I_normbyH(isnan(I_normbyH))=0;
I_inst_normbyH(isnan(I_inst_normbyH))=0;

%matrix of dominant links
for Target =1:ntargs
    for Source = 1:nvars
        
        vect_norm = reshape(I_normbyH(:,Source,Target),1,nlags); %normalized I (bits/bit)
        vect_nonnorm = reshape(I_lags(:,Source,Target),1,nlags);   %non-normalized I (bits)
        
        if Source==Target && mi.ZeroLagOpt ==1
            vect_norm(1)=0; %don't allow for zero-lag self links
            vect_nonnorm(1)=0;
        end
        
        if Source==Target && mi.NoSelfOpt==1 %skip self links if NoSelfOpt=1
            continue
        end

        if max(vect_norm)>0 && DomNormOpt == 1 %define dominant links based on normalized values
            [I_dom_normbyH(Source,Target), ind_lag] = max(vect_norm);    %maximum strength
            I_tau(Source,Target) = lagvect(ind_lag(1));
            TE_T(Source,Target) = TE(ind_lag(1),Source,Target);
            TE_T_normbyItot(Source,Target) = TE_normbyItot(ind_lag(1),Source,Target);

            
            I_dom(Source,Target) = I_lags(ind_lag,Source,Target);
        elseif max(vect_nonnorm>0) %define dominant links based on I in bits (non-normalized)
            [I_dom(Source,Target), ind_lag] = max(vect_nonnorm);    %maximum strength
            I_tau(Source,Target) = lagvect(ind_lag(1));
            TE_T(Source,Target) = TE(ind_lag(1),Source,Target);
            TE_T_normbyItot(Source,Target) = TE_normbyItot(ind_lag(1),Source,Target);
            
            I_dom_normbyH(Source,Target) = I_normbyH(ind_lag,Source,Target);
        end            
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
entropy.TE_tau=             TE;                 %transfer entropy TE(lag, source, target)
entropy.TE  =               TE_T;               %transfer entropy asociated with I_dom (TE at lag of strongest link)
entropy.TE_normbyItot =     TE_T_normbyItot;    %transfer entropy normalized by Itot = I(Source;Target_history;Target)

entropy.N_sources =         N_sources;          %number of sources detected to each target (same node can transmit at mult lags)
entropy.N_sourcenodes =     N_sourcenodes;      %number of individual source nodes to each target
entropy.I_inst =            I_inst;             %lag zero mutual information
entropy.I_inst_normbyH=     I_inst_normbyH;     %normalized lag zero mutual information
entropy.I_inst_sigthresh =  I_inst_sig;   %statistical significance threshold for mutual information

entropy.TE_Bonneville =     TE_Bonneville;  %AEG added for Nick salmon thesis

end

