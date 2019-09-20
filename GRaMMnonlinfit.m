function []=GRaMMnonlinfit(dataA,dataB)
%% GRaMMnonlinfit is used to builds the fitting curves of 
%%  nonlinear correlation pairs based on the common nonlinear models.
% Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM), 
%   for inter-correlations discovery among metabolome and microbiome. 
%   GRaMM considers the different characteristics of omics data, the effect of covariates, 
%   and the balance of linear and nonlinear correlations, by integrating the classical 
%   linear regression, the emerging maximum information coefficient (MIC).
% Authors :Dandan Liang, Mengci Li, Wei Jia, Tianlu Chen.
% Time: 2019.06

%% input :
% dataA: The Variable matrix1,such as the preprocessing metabolome data matrix.
% dataB: The Variable matrix2,such as the preprocessing microbiome data matrix.
%% five models
%model1: lin2p: 
% y = b(1).* x.^2 + b(2).* x + b(3);
%modeln2: log2p:
% y = b(1) .* log(x) + b(2);
%model3: exp2p:
% y = exp(b(1)+b(2) .* x);
% model4£¬s-curve function;
% y = 1./(b(1) + b(2).*exp(-x))
% model 5,power function,
% y = b(1).*x.^b(2)+b(3)

 warning('off')   
for i=1:size(dataA,1);
for j=1:size(dataB,1);
y = dataA(i,:); 
x = dataB(j,:);
[m,~]=max(x);[n,~]=min(x);x1=n:0.1:m;

%function1: lin2p: 
% y = b(1).* x.^2 + b(2).* x + b(3);
lin2p = inline('b(3).* x.^2 + b(2).* x + b(1)','b','x');
k = polyfit(x,y,2);
[beta1,r1,~]=nlinfit(x',y',lin2p ,k);
R1 = sum(r1.^2);
yhat1 = beta1(3).* x1.^2 + beta1(2).* x1 + beta1(1);

%function2: log2p:
 % y = b(1) .* log(x) + b(2);
 try
 syms l k 
equ1 = l.*log(x(1))+k -y(1) == 0;
equ2 = l.*log(x(2))+k -y(2) == 0;
[l,k] = solve(equ1,equ2,l,k);
l = vpa(l,2);
k = vpa(k,2);
syml = [l,k];
b2 = double(syml);
log2p = inline('b(1) .* log(x) + b(2)','b','x');
[beta2,r2,~]=nlinfit(x',y',log2p ,b2');
R2 = sum(r2.^2);
 yhat2 = beta2(1) .* log(x1) + beta2(2);
 catch
     R2=[];
 end

%funciton3: exp2p:
% y = exp(b(1)+b(2) .* x);
try
syms z f 
equ1 = log(y(1))-z-f.*x(1) == 0;
equ2 = log(y(2))-z-f.*x(2) == 0;
[z,f] = solve(equ1,equ2,z,f);
z = vpa(z,2);
f = vpa(f,2);
syme = [z,f];
b3 = double(syme);
exp2p = inline('exp(b(1)+b(2) .* x)','b','x');
[beta3,r3,~]=nlinfit(x',y',exp2p ,b3');
R3 = sum(r3.^2);
yhat3 = exp(beta3(1)+beta3(2) .* x1);
catch
    R3=[];
end

 
 % model5£¬s-curve function;
 %y = 1./(a + b.*exp(-x))
 try
syms u v 
equ1 = y(1)-1./(u + v.*exp(-x(1))) == 0;
equ2 = y(2)-1./(u + v.*exp(-x(2))) == 0;
[u, v] = solve(equ1,equ2,u,v);
u = vpa(u,2);
v = vpa(v,2);
syml = [u,v];
b5 = double(syml);
% curve-fitting;
s_fun5 = inline ('1./(b(1)+b(2).*exp(-x))','b','x');
% b5 = [1,1];
[beta5,r5,~]= nlinfit(x',y',s_fun5 ,b5);
R5 = sum(r5.^2);
yhat5 = 1./(beta5(1)+beta5(2).*exp(-x1));
 catch
     R5=[];
 end
 
 % model 6,power function,
 % y = b(1).*x.^b(2)+b(3)
 try
 syms w g r
equ1= w * x(1).^g +r - y(1) == 0;
equ2= w * x(2).^g +r - y(2) == 0;
equ3= w * x(3).^g +r - y(3) == 0;
[w,g,r] = solve(equ1,equ2,equ3,w,g,r);
w = vpa(w,2);
g = vpa(g,2);
r = vpa(r,2);
sympp = [w,g,r];
b6 = double(sympp);
power_fun6 = inline ('b(1).*x.^b(2)+b(3)','b','x');
[beta6,r6,~]= nlinfit(x',y',power_fun6 ,b6);
R6 = sum(r6.^2);
yhat6 = beta6(1).*x1.^beta6(2)+beta6(3);
 catch
     R6=[];
 end

 % graphic
 str=[ 'firM-' num2str(i) '* '  'sedM-' num2str(j)];
 filename=[ 'firM-' num2str(i)  'sedM-' num2str(j)];
 figure('name',str)
 [~,ii]=min([R1,R2,R3,R5,R6]) ;
 if ii==1
    plot(x,y,'.',x1,yhat1,'r')
    title (['y=' num2str(beta1(3)) '*x.^2+' num2str(beta1(2)) '*x+' num2str(beta1(1)) ';  SSE=' num2str(R1) ])
    set(gca,'FontName','Times New Roman','FontSize',8)
 elseif ii==2
    plot(x,y,'.',x1,yhat2,'r')
    title(['y=' num2str(beta2(1)) '*log(x)+' num2str(beta2(2)) '; SSE=' num2str(R2)])
    set(gca,'FontName','Times New Roman','FontSize',8)
 elseif ii==3
     plot(x,y,'.',x1,yhat3,'r')
    title(['y=exp(' num2str(beta3(1)) '+' num2str(beta3(2)) '.*x )' '; SSE=' num2str(R3)])
    set(gca,'FontName','Times New Roman','FontSize',8)
 elseif ii==4
    plot(x,y,'.',x1,yhat5,'r')
    title(['y=1./(' num2str(beta5(1)) '+' num2str(beta5(2)) '.*exp(-x) )' '; SSE=' num2str(R5)])
    set(gca,'FontName','Times New Roman','FontSize',8)  
 elseif ii==5
    plot(x,y,'.',x1,yhat6,'r');
    title(['y=' num2str(beta6(1)) '.*x.^' num2str(beta6(2)) '+' num2str(beta6(3)) '; SSE=' num2str(R5)])
    set(gca,'FontName','Times New Roman','FontSize',8) 
 end
 saveas(gcf,['D:\MATLAB.code\figure',filename,'.jpg']);       
end
end
 
 
 
 
 
   