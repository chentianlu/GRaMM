function [result,GRaMM] = GRaMM(A,B,C,Anorm,Bnorm,Brarefy,alpha,linR)
%% Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM), 
%  for inter-correlations discovery among metabolome and microbiome. 
%  GRaMM considers the different characteristics of omics data, the effect of covariates, 
%  and the balance of linear and nonlinear correlations, by integrating the classical 
%  linear regression, the emerging maximum information coefficient (MIC).
% Authors :Dandan Liang, Mengci Li, Wei Jia, Tianlu Chen.
% Time:2019.06
    %% input:
%      For the inputted matrix, raws are variates, columns are samples.
% A is a matrix of metabolome data, p*n. The raw(p) is the variates,column(n) is samples.
% B is a matrix of microbiome data,v*n.The raw(v) is the variates,column(n) is samples.
% C is optional, presents the covariates matrix£¬u*n. The raw(u) is the covariates,column(n) is samples.
% Anorm = 'yes',or 'no'(defualt). If choose "yes", the metabolome data will
%       be benormalized by samples."no" means the metabolome data will not be
%       normalized by samples.
% Bnorm = 'yes',or 'no'(defualt). If choose "yes", the microbiome data will
%       be benormalized by samples."no" means the microbiome data will not be
%       normalized by samples.
% Brarefy = 'yes',or 'no'(defualt). If choose "yes", the microbiome data will
%        be normalized using rarefying.
% alpha(defualt = 0.05),the threshold of significant difference for linear correlation.
% linR (defualt = 0.1). The linear correlation threshold to determine whether or not it is linear.
%        "> linR" means it tends to a strong linear correlation;
%       " <linR " means it it's more llikely to be nonlinear correlation.
    %% output
% "results1" includes 4 columns, correlation strength(r),p value,FDR value and
%       the correlation types("0"=linear correlation, "1"=nonlinear correlation).
% "results2" includes 4 matrixes(v*p) for r,p,FDR and type respectively.

warning('off') 
if nargin <= 7,linR= 0.1;end
if nargin <= 6,alpha = 0.05;end
if nargin <= 5,Brarefy='no';end
if nargin <= 4, Bnorm = 'no';end
if nargin <= 3,Anorm = 'no';end
if nargin<=2,C=[];end
 
%% data preprocessing analysis
% NAN and zero filling for metabolomic data using KNN, defualt k=10;
A(A==0)=nan;
A=knnimpute(A, 10);
%rarefying for microbial data;
switch Brarefy
    case 'yes'
       B =rarefying(B);
    case 'no'
       B = B;
end

% count normaliztion for metabolomic data and microbiome data;
switch Anorm
case 'yes'
[~,nn] = size(A);
dA = [];
for kk = 1:nn;
    sumA = sum(A(:,kk));
    norA = 30000.* A(:,kk)/sumA;
    dA = [dA,norA];
end
dataA = log( dA+ 1);
case 'no'
    dataA = log( A+ 1);
end

switch Bnorm
case 'yes'
    [m,n] = size(B);
dataB = [];
dB = [];
for kn = 1:n;
    sumB = sum(B(:,kn));
    norB = 30000.* B(:,kn)/sumB;
    dB = [dB,norB];
end
index = dB == 0;
dB(index)=1;
for k = 1:n;
    gmB = (prod(dB(:,k))^(1/m));  
    logB = log(dB(:,k))/log(gmB);  
    dataB = [dataB,logB];
end
case 'no'
[m,n] = size(B);
dataB = [];
index = B == 0;
B(index)=1;
for k = 1:n;
    gmB = (prod(B(:,k))^(1/m));  
    logB = log(B(:,k))/log(gmB);  
    dataB = [dataB,logB];
end
end
%--------------------------------------------------------------------------
%% correlation analysis:
result =[];
for i=1:size(dataA,1);
for j=1:size(dataB,1);
y = dataA(i,:); 
x1 = dataB(j,:);
x2 = C; 
X = [ones(size(x1));x1;x2];
[b,~,~,~,~] = regress(y',X');
beta1 = b(2); 
r1 = beta1*(std(x1)/std(y)); %Calculate the r value corresponding to the partial regression coefficient.
 [~,pvalue,~,~,~] = regresslm([x1;x2]',y');
 p = pvalue(2);
 absr=abs(r1);

if  p < alpha && absr > linR
    regress_re = [r1,p,0];
    result = [result;regress_re];
   
else
    empty = isempty(C); perm_n=100;
    if empty == 0; 
       beta2 = b(3); 
       m_y = y-(beta2*x2);
       y = dataA(i,:); x1 = dataB(j,:); x2 = C; 
       Y = [ones(size(y));y; x2];  
       [xb,~,~,~,~] = regress(x1',Y'); 
       d2 = xb(3);  
       m_x = x1-(d2*x2);
       
       statsmic = mine (m_x,m_y);
       mic = statsmic.mic;
       mun_p =0;
    for k = 1:100;
        perm_y = m_y(randperm(length(y)));
        perm_statsmic = mine (perm_y,m_x);
        perm_mic = perm_statsmic.mic;
        if mic <= perm_mic;
           mun_p = mun_p+1 ;
        end
        mic_p = mun_p/(perm_n+1);
    end
    else 
       statsmic = mine (x1,y); 
       mic = statsmic.mic; 
       mun_p = 0;
    for  k = 1:100 
        perm_y = y(randperm(length(y)));
        perm_statsmic = mine (x1,perm_y);
        perm_mic = perm_statsmic.mic;
        if mic <= perm_mic;
           mun_p = mun_p + 1 ;
        end
        mic_p = (mun_p+1) ./ (perm_n + 1);
    end
    end
   %----------------------------------------------------------------- 
    mic_re = [mic,mic_p,1];
    result =[result;mic_re];
end 
end
end
rawP = result(:,2);
FDR = mafdr(rawP,'BHFDR', true);
result=[result,FDR];
GRaMM.FDR = reshape(FDR,j,i);
GRaMM.r = reshape(result(:,1),j,i);
GRaMM.p = reshape(result(:,2),j,i);
GRaMM.type = reshape(result(:,4),j,i);
end

%--------------------------------------------------------------------------
function[b,p,F,t,R2]=regresslm(x,y)
%Calculate the normalized values of r and p.
warning('off');
[n,k]=size(x);
X=[ones(n,1),x];
A=X'*X; 
C=inv(A); 
b=X\y; 
RSS=y'*y-b'*X'*y; 
MSe=RSS/(n-k-1);
Up=b.*b./diag(C);
F=Up/MSe;
sb=sqrt(MSe*diag(C)); 
t=b./sb; 
SSy=var(y)*(n-1);
R2=(SSy-RSS)/SSy;
p = 2*(1-tcdf(abs(t),n-k-1));
end


function[rarefy_dataB]=rarefying(dataB)

dataB=round(dataB);
sumc=sum(dataB); 
min_col=ceil(min(min(sumc)));
[a,b]=size(dataB);
re_dataB=[];
for j=1:b
new_otu=[];
for i= 1:a
otui_allname =  i; counti=dataB(i,j);
repotui= repelem(otui_allname,counti );
new_otu=[new_otu,repotui];
end
re_otu=randsample(new_otu,min_col);
rep_otu=[];
for i=1:a
    new_counti=sum(re_otu==i);
    rep_otu=[rep_otu,new_counti];
    
end
re_dataB=[re_dataB;rep_otu];
end
rarefy_dataB = re_dataB';
end




    

   

    
    
    
    
    
    
    
    
    
    
    
    
  
