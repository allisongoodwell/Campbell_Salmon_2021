function [pdf, Edges]= compute_pdfGUI(Data,N,bin_scheme,Range)
% Compute pdf 1D, 2D, or 3D
% fixed binning only

%N = number of grid points at which to compute kernel
%bin_scheme: local or global
%Range: used for global binning scheme only

dim = size(Data,2);

if dim ==1
Data(isnan(Data))=[]; 
else
nancheck = sum(isnan(Data')); 
%disp(nancheck)
%if a tuple has a nan value, ignore it
Data(nancheck>0,:)=[];
end

%disp(size(Data))

nTup = size(Data,1);


if dim==1
    pdf = zeros(1,N);
elseif dim==2
    pdf=zeros(N,N);
elseif dim==3
    pdf=zeros(N,N,N);
end

Coords = zeros(dim,N);
Edges = zeros(dim,N+1);
xo = zeros(1,dim);

%%%%%%%%%%%%%%%%%%%%% Determine bin coordinates %%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:dim
    
    if strcmp(bin_scheme,'local')==1
        Edges(i,:)= linspace(min(Data(:,i)-10^-8), max(Data(:,i)),N+1);
        xo(i) = min(Data(:,i));
    elseif strcmp(bin_scheme,'global')==1
        Edges(i,:) = linspace(Range(1,i)-10^-8,Range(2,i),N+1);
        xo(i) = Range(1,i);
    end
    Coords(i,:)=(Edges(i,1:end-1)+Edges(i,2:end))./2;
    %Coords(i,1)=0;
    
end


for i = 1:dim
    dat = Data(:,i);
    bindata = ones(size(dat));
    edges = Edges(i,:);

    for e = 1:N
            bindata(dat >= edges(e) & dat<edges(e+1))= e;
    end
    
    BinData(:,i)=bindata; %(i,j,k) bin numbers for each data point
end



if dim==1
    C=zeros(1,N);
    for n = 1:nTup
        dat = BinData(n);
        C(dat)=C(dat)+1;
    end
    
elseif dim==2
    C=zeros(N,N);
    for n = 1:nTup
        dat = BinData(n,:);
        C(dat(1),dat(2))=C(dat(1),dat(2))+1;
    end
elseif dim==3
    C=zeros(N,N,N);
    for n = 1:nTup
        dat = BinData(n,:);
        C(dat(1),dat(2),dat(3))=C(dat(1),dat(2),dat(3))+1;
    end
end

pdf = C./nTup;

end



