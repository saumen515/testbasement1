 % % *** and SVM model                                                                ***
% % *** To run the script Netlab and LSSVM softwares/toolboxes are required to be 
% % *** downloaded in the same folder. 
% % *** Tested using Netlab software developed by I. Nabney, [Nabney, I., 2002. NETLAB: 
% % *** algorithms for pattern recognition. Springer Science & Business
% % *** Media]; Netlab can be downloaded at https://www2.aston.ac.uk/eas/research/groups/ncrg/resources/netlab/downloads
% % *** Tested using LSSVM software/toolbox which cab be downladed at http://www.esat.kuleuven.be/sista/lssvmlab/
% % *** Tested using Matlab R2017a under Windows 7 
% SVM Regression.......................................................
clc;
clear all;
close all;
 
%Training data;
load gravdeep.dat;
aa=gravdeep(:,1);%longitude
bb=gravdeep(:,2);%latitude
cc=gravdeep(:,3);%BG value
 
dd=gravdeep(:,4);%Altitude
ee=gravdeep(:,5);%deep basement
% %----------------------------------------------------
 
totalden = aa;
totalporo = bb;
totalgammar = cc;
%........................................
 
%Noise stability test
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
Aa1=0.99;%0.99
dsm(1)=a(1);
for i=2:nmax,
dsm(i)=Aa1.*dsm(i-1)+r(i);
end
dsm=dsm/max(abs(dsm));
dsm=dsm';
%r=r';
%den0=den.*r;
%den=a;
%totalden0=dsm;
totalden0=totalden.*dsm;
totalden10=totalden+0.001.*totalden0;
totalden20=totalden+0.002.*totalden0;
totalden30=totalden+0.003.*totalden0;
totalden40=totalden+0.004.*totalden0;
totalden50=totalden+0.005.*totalden0;
 
randn('state',0);
nmax=211;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
pr(1)=r(1);
Ba2=-0.0957;%-0.0957
for i=2:nmax,
pr(i)=Ba2.*pr(i-1)+r(i);
end
 
pr=pr/max(abs(pr));
pr=pr';
%poro0=poro.*r;
%poro=b;
%totalporo0=pr;
totalporo0=totalporo.*pr;
totalporo10=totalporo+0.001.*totalporo0;
totalporo20=totalporo+0.002.*totalporo0;
totalporo30=totalporo+0.003.*totalporo0;
totalporo40=totalporo+0.004.*totalporo0;
totalporo50=totalporo+0.005.*totalporo0;
 
randn('state',0);
nmax=211;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
gg(1)=r(1);
Ca3=0.0142;%0.0142
for i=2:nmax,
gg(i)=Ca3.*gg(i-1)+r(i);
end
gg=gg/max(abs(gg));
gg=gg';
%gammar0=gammar.*r;
%gammar=gammar';
%totalgammar0=gg;
totalgammar0=totalgammar.*gg;
totalgammar10=totalgammar+0.001.*totalgammar0;
totalgammar20=totalgammar+0.002.*totalgammar0;
totalgammar30=totalgammar+0.003.*totalgammar0;
totalgammar40=totalgammar+0.004.*totalgammar0;
totalgammar50=totalgammar+0.005.*totalgammar0;
 
fid = fopen('totalinerror10.dat','w');
totalerror10=[totalden10';totalporo10';totalgammar10'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror10);
fclose(fid);
 
fid = fopen('totalinerror20.dat','w');
totalerror20=[totalden20';totalporo20';totalgammar20'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror20);
fclose(fid);
 
fid = fopen('totalinerror30.dat','w');
totalerror30=[totalden30';totalporo30';totalgammar30'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror30);
fclose(fid);
 
fid = fopen('totalinerror40.dat','w');
totalerror40=[totalden40';totalporo40';totalgammar40'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror40);
fclose(fid);
 
fid = fopen('totalinerror50.dat','w');
totalerror50=[totalden50';totalporo50';totalgammar50'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror50);
fclose(fid);
 
load totalinerror10.dat
load totalinerror20.dat
load totalinerror30.dat
load totalinerror40.dat
load totalinerror50.dat
 
%Training data;
totalinerror=totalinerror30;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
ccn=totalinerror(:,3);%BG value
 
%........................................................................
% Real data
load gravwhole.dat;
ff=gravwhole(:,1);%longitude
gg=gravwhole(:,2);%latitude
hh=gravwhole(:,3);%BG value
ii=gravwhole(:,4);%altitude
 
x=[aan,bbn];
t=ee;
 
x=x';
t=t';
 
%.....................................................noise adding......
 %ndata=211;
 noise = 0.15;          % Standard deviation of noise distribution.
 
 %*********************************************************
target=t;
randn('state',0);
nmax=211;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
Ta1=0.99;%0.99
dm(1)=a(1);
for i=2:nmax,
dm(i)=Ta1.*dm(i-1)+r(i);
end
 
dm=dm/max(abs(dm));
target0=target.*dm;
target10=target+0.1.*target0;
target20=target+0.2.*target0;
target30=target+0.3.*target0;
target40=target+0.4.*target0;
target50=target+0.5.*target0;
t=target30;
 
%**********************************************************************
 
%***********************************************************************
[xn,minp,maxp,tn,mint,maxt] = premnmx(x,t);
[xtrans,transMat] = prepca(xn,0.0000000000000002);
[R,Q] = size(xtrans);
iitst = 2:4:Q;
iival = 4:4:Q;
iitr = [1:4:Q 3:4:Q];
valP = xtrans(:,iival); val.T = tn(:,iival);
testP = xtrans(:,iitst); test.T = tn(:,iitst);
xtr = xtrans(:,iitr); ttr = tn(:,iitr);
xr=xtr';
zr=ttr';
 
x=xr;
y=zr;
fprintf('Program paused. Press enter to continue.\n');
%pause;
 
%----------------Parameter specifications--------------------------------------
C = 120;
lambda = 0.01;
epsilon = 0.01;
kerneloption = 01;
 %kernel='poly';
kernel='Gaussian';
%kernel='htrbf';
 %kernel='wavelet';
 %kernel='polyhomog';
% kernel='frame';       
 
verbose=1;
[xsup,ysup,w,w0] = svmreg(x,y,C,epsilon,kernel,kerneloption,lambda,verbose);
%--------------------------------------------------------
 
 
xtest=xr;
[ypredtr] = svmval(xtest,xsup,w,w0,kernel,kerneloption);
[ypred] = postmnmx(ypredtr',mint,maxt);
ee=ee';
tr=ee(:,iitr);
figure
%,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
%error bar calculation with SVR
type = 'function approximation';
[Yp,alpha,b,gam,sig2] = lssvm(x,y,type);
sig2e = bay_errorbar({x,y,type,gam,sig2},'figure');
 
figure
xis=1:106;
E1=sig2(:,1).*ones(size(Yp));
%subplot(3,1,1),
errorbar(xis,ypred,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tr,'ro-',xis, ypred,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -6 6]);
legend('Error-bar','Target','SVMreg.');
xlabel('No. of Data','Fontsize',20);
ylabel('Sediment Depth(Km)','Fontsize',20);
title('Training Interval')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%_____________________
 
 
%............................................
%sqvalx=sigval;
valX=valP;
xtest=valX';
ypredtr = svmval(xtest,xsup,w,w0,kernel,kerneloption);
[ysvrval] = postmnmx(ypredtr',mint,maxt);
av=ysvrval;
tav=ee(:,iival);
%error bar calculation with SVR:Validation Interval
type = 'function approximation';
[Yp,alpha,b,gam,sig2] = lssvm(xtest,tav',type);
sig2e = bay_errorbar({x,y,type,gam,sig2},'figure');
 
figure
xis=1:52;
E1=sig2(:,1).*ones(size(Yp));
%subplot(3,1,1),
errorbar(xis,ysvrval,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tav,'ro-',xis, ysvrval,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -6 6]);
legend('Error-bar','Target','SVMreg.');
xlabel('No. of Data','Fontsize',20);
ylabel('Sediment Depth(Km)','Fontsize',20);
title('Validation Interval')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%sqtestx=sigtst;
testX=testP;
xtest=testX';
ypredtr = svmval(xtest,xsup,w,w0,kernel,kerneloption);
[ysvrtst] = postmnmx(ypredtr',mint,maxt);
ate=ysvrtst;
% ate = mlpfwd(net,testP');
% [ate] = postmnmx(ate',mint,maxt);
tat=ee(:,iitst);
%error bar calculation with SVR:Validation Interval
type = 'function approximation';
[Yp,alpha,b,gam,sig2] = lssvm(xtest,tat',type);
sig2e = bay_errorbar({x,y,type,gam,sig2},'figure');
figure
xis=1:53;
E1=sig2(:,1).*ones(size(Yp));
%subplot(3,1,1),
errorbar(xis,ysvrtst,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tat,'ro-',xis, ysvrtst,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -6 6]);
legend('Error-bar','Target','SVMreg.');
xlabel('No. of Data','Fontsize',24);
ylabel('Sediment Depth(Km)','Fontsize',24);
title('Test Interval')
set(gca,'LineWidth',5,'FontName','Tahoma','Fontsize',24);axis tight
% figure(28);
;
% 
figure
% Real data analysis;
vinerror10=[ff,gg];
[e10]= tramnmx(vinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
a10 = svmval(p2e10',xsup,w,w0,kernel,kerneloption);
%a10 = gpfwd(net,p2e10');
[a10] = postmnmx(a10',mint,maxt);
%[a10]=-a10;
 
fdata=488;
gdata=488;
fmax=max(ff);
gmax=max(gg);
fmin=min(ff);
gmin=min(gg);
fv = linspace(fmin, fmax,fdata);
gv = linspace(gmin, gmax,gdata);
[fi,gi] = meshgrid(fv,gv);
c1i = griddata(ff,gg,a10,fi,gi,'cubic');
 
meshz(fi,gi,c1i);hold on;
colorbar;
colormap(jet);
%meshc(fi,gi,c1i);hold on;
%surf(fi,gi,c1i);hold on;
%surf(ai,bi,c1i);hold on;
%[cout,H,cf]=contourf(fi,gi,c1i);
%colorbarf(cout,H);hold on
plot3(ff,gg,a10,'v','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
%view(10,90);
% view(0,90);
%shading interp;
grid on;
box on;
%axis tight;
%axis auto;
%axis square
colorbar;
set(gca,'LineWidth',2,'FontName','Tahoma','Fontsize',20);axis tight;
xlabel('Longitude(in degree)','Fontsize',20);
ylabel('Latitude(in degree)','Fontsize',20);
zlabel('Deeper Interface Depth(Km)','Fontsize',20);set(gca,'Zdir','reverse');
title('Deeper Interface Depth Estimation','Fontsize',20);
%.........................................................................
 
 
Training....................................................
%.........................................................................
ntotal=211;
ntrain=106;
nval=52;
ntest=53;
 
x1train=tr(1:ntrain);
x1val=tav(1:nval);
x1test=tat(1:ntest);
 
x2train=ypred(1:ntrain);
x2val=av(1:nval);
x2test=ate(1:ntest);
 
disp('Statistics of x1/obs......................................');
x1avgtrain=mean(x1train)
x1avgnval=mean(x1val)
x1avgtest=mean(x1test)
 
x1stdtrain=std(x1train)
x1stdval=std(x1val)
x1stdtest=std(x1test)
 
x1trmedian=median(x1train)
x1valmedian=median(x1val)
x1testmedian=median(x1test)
 
x1trmode=mode(x1train)
x1valmode=mode(x1val)
x1testmode=mode(x1test)
 
x1trskewness=skewness(x1train)
x1valskewness=skewness(x1val)
x1testskewness=skewness(x1test)
 
x1trkurtosis=kurtosis(x1train)
x1valkurtosis=kurtosis(x1val)
x1testkurtosis=kurtosis(x1test)
 
disp('Statistics of x2/SVMreg......................................')
x2avgtrain=mean(x2train)
x2avgnval=mean(x2val)
x2avgtest=mean(x2test)
 
x2stdtrain=std(x2train)
x2stdval=std(x2val)
x2stdtest=std(x2test)
 
x2trmedian=median(x2train);
x2valmedian=median(x2val)
x2testmedian=median(x2test)
 
x2trmode=mode(x2train)
x2valmode=mode(x2val)
x2testmode=mode(x2test)
 
x2trskewness=skewness(x2train)
x2valskewness=skewness(x2val)
x2testskewness=skewness(x2test)
 
x2trkurtosis=kurtosis(x2train)
x2valkurtosis=kurtosis(x2val)
x2testkurtosis=kurtosis(x2test)
 
x2errtrain=mse(x2train-x1train)
x2errval=mse(x2val-x1val)
x2errtest=mse(x2test-x1test)
 
x2maetrain=mae(x2train-x1train)
x2maeval=mae(x2val-x1val)
x2maetest=mae(x2test-x1test)
 
x2retrain=1.0-sum((x1train-x2train).^2)./sum(x1train.^2)
x2reval=1.0-sum((x1val-x2val).^2)./sum(x1val.^2)
x2retest=1.0-sum((x1test-x2test).^2)./sum(x1test.^2)
 
x2dtrain=1.0-sum((x1train-x2train).^2)./sum(abs(x2train-mean(x1train))+abs(x1train-mean(x1train)).^2)
x2dval=1.0-sum((x1val-x2val).^2)./sum(abs(x2val-mean(x1val))+abs(x1val-mean(x1val)).^2)
x2dtest=1.0-sum((x1test-x2test).^2)./sum(abs(x2test-mean(x1test))+abs(x1test-mean(x1test)).^2)
 
trcorrcoef=corrcoef(x1train,x2train)
valcorrcoef=corrcoef(x1val,x2val)
testcorrcoef=corrcoef(x1test,x2test)
 
trcod=diag(trcorrcoef,1).*diag(trcorrcoef,1)
valcod=diag(valcorrcoef,1).*diag(valcorrcoef,1)
testcod=diag(testcorrcoef,1).*diag(testcorrcoef,1)
 
%***************Training Error Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error analysis**********************************************************************
aa=aa';
bb=bb';
cc=cc';
 
inputden = aa(:,iitr);
inputporo = bb(:,iitr);
inputgammar = cc(:,iitr);
 
inputden=inputden';
inputporo=inputporo';
inputgammar=inputgammar';
 
randn('state',0);
nmax=106;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
a1=0.94;%0.94
dtr(1)=a(1);
for i=2:nmax,
dtr(i)=a1.*dtr(i-1)+r(i);
end
 
dtr=dtr/max(abs(dtr));
dtr=dtr';
%r=r';
%den0=den.*r;
%den=a;
%inputden0=dtr;
inputden0=inputden.*dtr;
inputden10=inputden+0.001.*inputden0;
inputden20=inputden+0.002.*inputden0;
inputden30=inputden+0.003.*inputden0;
inputden40=inputden+0.004.*inputden0;
inputden50=inputden+0.005.*inputden0;
 
randn('state',0);
nmax=106;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
ptr(1)=r(1);
a2=-0.107;%-0.107
for i=2:nmax,
ptr(i)=a2.*ptr(i-1)+r(i);
end
ptr=ptr';
ptr=ptr/max(abs(ptr));
%poro0=poro.*r;
%poro=b;
%inputporo0=ptr;
inputporo0=inputporo.*ptr;
inputporo10=inputporo+0.001.*inputporo0;
inputporo20=inputporo+0.002.*inputporo0;
inputporo30=inputporo+0.003.*inputporo0;
inputporo40=inputporo+0.004.*inputporo0;
inputporo50=inputporo+0.005.*inputporo0;
 
randn('state',0);
nmax=106;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
ggtr(1)=r(1);
a3=0.122;%0.122
for i=2:nmax,
ggtr(i)=a3.*ggtr(i-1)+r(i);
end
 
ggtr=ggtr/max(abs(ggtr));
ggtr=ggtr';
%gammar0=gammar.*r;
%gammar=gammar';
%inputgammar0=ggtr;
inputgammar0=inputgammar.*ggtr;
inputgammar10=inputgammar+0.001.*inputgammar0;
inputgammar20=inputgammar+0.002.*inputgammar0;
inputgammar30=inputgammar+0.003.*inputgammar0;
inputgammar40=inputgammar+0.004.*inputgammar0;
inputgammar50=inputgammar+0.005.*inputgammar0;
 
fid = fopen('sinerror10.dat','w');
%inputerror10=[inputden10';inputporo10';inputgammar10'];
inputerror10=[inputden10';inputporo10'];
fprintf(fid,'%6.4f %6.4f\n',inputerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror10);
fclose(fid);
 
fid = fopen('sinerror20.dat','w');
%inputerror20=[inputden20';inputporo20';inputgammar20'];
inputerror20=[inputden20';inputporo20'];
fprintf(fid,'%6.4f %6.4f\n',inputerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror20);
fclose(fid);
 
fid = fopen('sinerror30.dat','w');
%inputerror30=[inputden30';inputporo30';inputgammar30'];
inputerror30=[inputden30';inputporo30'];
fprintf(fid,'%6.4f %6.4f\n',inputerror30);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror30);
fclose(fid);
 
fid = fopen('sinerror40.dat','w');
%inputerror40=[inputden40';inputporo40';inputgammar40'];
inputerror40=[inputden40';inputporo40'];
fprintf(fid,'%6.4f %6.4f\n',inputerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror40);
fclose(fid);
 
fid = fopen('sinerror50.dat','w');
%inputerror50=[inputden50';inputporo50';inputgammar50'];
inputerror50=[inputden50';inputporo50'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror50);
fprintf(fid,'%6.4f %6.4f\n',inputerror50);
fclose(fid);
 
load sinerror10.dat;
load sinerror20.dat;
load sinerror30.dat;
load sinerror40.dat;
load sinerror50.dat
 
 
[e10tr]= tramnmx(sinerror10',minp,maxp);
[p2e10tr] = trapca(e10tr,transMat);
a10tr = svmval(p2e10tr',xsup,w,w0,kernel,kerneloption);
%a10 = gpfwd(net,p2e10');
[a10tr] = postmnmx(a10tr',mint,maxt);
 
 
 
%t=t(:,iitr);
%tr=t(:,iitr);
trdev10=tr-a10tr;
trdev10=trdev10';
trdp10=trdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
trmse10=mse(trdp10);
 
figure(29);
 
subplot(1,1,1),plot(trdp10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e20tr]= tramnmx(sinerror20',minp,maxp);
[p2e20tr] = trapca(e20tr,transMat);
a20tr = svmval(p2e20tr',xsup,w,w0,kernel,kerneloption);
%a20tr = gpfwd(net,p2e20tr');
[a20tr] = postmnmx(a20tr',mint,maxt);
 
 
 
trdev20=tr-a20tr;
trdev20=trdev20';
trdp20=trdev20(:,1);
trmse20=mse(trdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(30);
 
subplot(1,1,1),plot(trdp20,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e30tr]= tramnmx(sinerror30',minp,maxp);
[p2e30tr] = trapca(e30tr,transMat);
a30tr = svmval(p2e30tr',xsup,w,w0,kernel,kerneloption);
%a30tr = gpfwd(net,p2e30tr');
[a30tr] = postmnmx(a30tr',mint,maxt);
 
 
 
trdev30=tr-a30tr;
trdev30=trdev30';
trdp30=trdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
trmse30=mse(trdp30);
 
figure(31);
 
subplot(1,1,1),plot(trdp30,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e40tr]= tramnmx(sinerror40',minp,maxp);
[p2e40tr] = trapca(e40tr,transMat);
a40tr = svmval(p2e40tr',xsup,w,w0,kernel,kerneloption);
%a40tr = gpfwd(net,p2e40tr');
[a40tr] = postmnmx(a40tr',mint,maxt);
 
trdev40=tr-a40tr;
trdev40=trdev40';
trdp40=trdev40(:,1);
trmse40=mse(trdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(32);
 
subplot(1,1,1),plot(trdp10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e50tr]= tramnmx(sinerror50',minp,maxp);
[p2e50tr] = trapca(e50tr,transMat);
a50tr = svmval(p2e50tr',xsup,w,w0,kernel,kerneloption);
%a50tr = gpfwd(net,p2e50tr');
[a50tr] = postmnmx(a50tr',mint,maxt);
 
 
 
trdev50=tr-a50tr;
trdev50=trdev50';
trdp50=trdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
trmse50=mse(trdp50);
 
figure(33);
subplot(1,1,1),plot(trdp50,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%error diagram
figure(34)
 
nx=[10 20 30 40 50];
msetr=[trmse10 trmse20 trmse30 trmse40 trmse50];
bar(nx,msetr);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Training Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Obs.','SVMreg', 'SVMreg+STD','SVMreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%******************Validation Error ANALYSIS****************************
 
 
vdenv = aa(:,iival);
vporov = bb(:,iival);
vgammarv = cc(:,iival);
 
valden=vdenv';
valporo=vporov';
valgammar=vgammarv';
 
randn('state',0);
nmax=52;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
dv(1)=r(1);
A1=0.99;%0.99
for i=2:nmax,
dv(i)=A1.*dv(i-1)+r(i);
end
 
dv=dv/max(abs(dv));
dv=dv';
%r=r';
%den0=den.*r;
%valden0=dv;
valden0=valden.*dv;
valden10=valden+0.001.*valden0;
valden20=valden+0.002.*valden0;
valden30=valden+0.003.*valden0;
valden40=valden+0.004.*valden0;
valden50=valden+0.005.*valden0;
 
randn('state',0);
nmax=52;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
pv(1)=r(1);
A2=-0.0495;%-0.0495
for i=2:nmax,
pv(i)=A2.*pv(i-1)+r(i);
end
pv=pv';
pv=pv/max(abs(pv));
%poro0=poro.*r;
%poro=b;
%valporo0=pv;
valporo0=valporo.*pv;
valporo10=valporo+0.001.*valporo0;
valporo20=valporo+0.002.*valporo0;
valporo30=valporo+0.003.*valporo0;
valporo40=valporo+0.004.*valporo0;
valporo50=valporo+0.005.*valporo0;
 
randn('state',0);
nmax=52;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
gva(1)=r(1);
A3=0.1408;%0.1408
for i=2:nmax,
gva(i)=A3.*gva(i-1)+r(i);
end
 
gva=gva/max(abs(gva));
gva=gva';
%gammar0=gammar.*r;
%gammar=gammar';
%valgammar0=gva;
valgammar0=valgammar.*gva;
valgammar10=valgammar+0.001.*valgammar0;
valgammar20=valgammar+0.002.*valgammar0;
valgammar30=valgammar+0.003.*valgammar0;
valgammar40=valgammar+0.004.*valgammar0;
valgammar50=valgammar+0.005.*valgammar0;
 
fid = fopen('vinerror10.dat','w');
%valerror10=[valden10';valporo10';valgammar10'];
valerror10=[valden10';valporo10'];
fprintf(fid,'%6.4f %6.4f\n',valerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror10);
fclose(fid);
 
fid = fopen('vinerror20.dat','w');
%valerror20=[valden20';valporo20';valgammar20'];
valerror20=[valden20';valporo20'];
fprintf(fid,'%6.4f %6.4f\n',valerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror20);
fclose(fid);
 
fid = fopen('vinerror30.dat','w');
%valerror30=[valden30';valporo30';valgammar30'];
valerror30=[valden30';valporo30'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror30);
fprintf(fid,'%6.4f %6.4f\n',valerror30);
fclose(fid);
 
fid = fopen('vinerror40.dat','w');
%valerror40=[valden40';valporo40';valgammar40'];
valerror40=[valden40';valporo40'];
fprintf(fid,'%6.4f %6.4f\n',valerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror40);
fclose(fid);
 
fid = fopen('vinerror50.dat','w');
%valerror50=[valden50';valporo50';valgammar50'];
valerror50=[valden50';valporo50'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror50);
fprintf(fid,'%6.4f %6.4f\n',valerror50);
fclose(fid);
 
load vinerror10.dat;
load vinerror20.dat;
load vinerror30.dat;
load vinerror40.dat;
load vinerror50.dat
 
 
[e10v]= tramnmx(vinerror10',minp,maxp);
[p2e10v] = trapca(e10v,transMat);
a10v = svmval(p2e10v',xsup,w,w0,kernel,kerneloption);
%a10v = gpfwd(net,p2e10v');
[a10v] = postmnmx(a10v',mint,maxt);
 
tav=tav';
vdev10=tav'-a10v;
vdev10=vdev10';
vdp10=vdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
valmse10=mse(vdp10);
 
figure(35);
 
subplot(1,1,1),plot(vdp10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e20v]= tramnmx(vinerror20',minp,maxp);
[p2e20v] = trapca(e20v,transMat);
a20v = svmval(p2e20v',xsup,w,w0,kernel,kerneloption);
%a20v = gpfwd(net,p2e20v');
[a20v] = postmnmx(a20v',mint,maxt);
 
vdev20=tav'-a20v;
vdev20=vdev20';
vdp20=vdev20(:,1);
valmse20=mse(vdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(36);
 
subplot(1,1,1),plot(vdp20,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e30v]= tramnmx(vinerror30',minp,maxp);
[p2e30v] = trapca(e30v,transMat);
a30v = svmval(p2e30v',xsup,w,w0,kernel,kerneloption);
%a30v = gpfwd(net,p2e30v');
[a30v] = postmnmx(a30v',mint,maxt);
 
 
 
vdev30=tav'-a30v;
vdev30=vdev30';
vdp30=vdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
valmse30=mse(vdp30);
 
figure(37);
 
subplot(1,1,1),plot(vdp30,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e40v]= tramnmx(vinerror40',minp,maxp);
[p2e40v] = trapca(e40v,transMat);
a40v = svmval(p2e40v',xsup,w,w0,kernel,kerneloption);
%a40v = gpfwd(net,p2e40v');
[a40v] = postmnmx(a40v',mint,maxt);
 
vdev40=tav'-a40v;
vdev40=vdev40';
vdp40=vdev40(:,1);
valmse40=mse(vdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(38);
 
subplot(1,1,1),plot(vdp40,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e50v]= tramnmx(vinerror50',minp,maxp);
[p2e50v] = trapca(e50v,transMat);
a50v = svmval(p2e50v',xsup,w,w0,kernel,kerneloption);
%a50v = gpfwd(net,p2e50v');
[a50v] = postmnmx(a50v',mint,maxt);
 
 
 
vdev50=tav'-a50v;
vdev50=vdev50';
vdp50=vdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
valmse50=mse(vdp50);
 
figure(40);
subplot(1,1,1),plot(vdp50,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%error diagram
figure(41)
nx=[10 20 30 40 50];
mseval=[valmse10 valmse20 valmse30 valmse40 valmse50];
bar(nx,mseval);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Validation Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Obs.','SVMreg', 'SVMreg+STD','SVMreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%**********************Test Error Analysis*****************************
%Error Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error analysis
dent = aa(:,iitst);
porot = bb(:,iitst);
gammart = cc(:,iitst);
 
testden=dent';
testporo=porot';
testgammar=gammart';
 
randn('state',0);
nmax=53;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
dtst(1)=r(1);
B1=0.99;%0.99
 
for i=2:nmax,
dtst(i)=B1.*dtst(i-1)+r(i);
end
 
dtst=dtst/max(abs(dtst));
dtst=dtst';
%r=r';
%den0=den.*r;
%den=a;
%testden0=dtst;
testden0=testden.*dtst;
testden10=testden+0.001.*testden0;
testden20=testden+0.002.*testden0;
testden30=testden+0.003.*testden0;
testden40=testden+0.004.*testden0;
testden50=testden+0.005.*testden0;
 
randn('state',0);
nmax=53;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
pt(1)=r(1);
B2=-0.2317;%-0.2317
for i=2:nmax,
pt(i)=B2.*pt(i-1)+r(i);
end
 
pt=pt/max(abs(pt));
pt=pt';
%poro0=poro.*r;
%poro=b;
%testporo0=pt;
testporo0=testporo.*pt;
testporo10=testporo+0.001.*testporo0;
testporo20=testporo+0.002.*testporo0;
testporo30=testporo+0.003.*testporo0;
testporo40=testporo+0.004.*testporo0;
testporo50=testporo+0.005.*testporo0;
 
randn('state',0);
nmax=53;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
gt(1)=r(1);
B3=-0.036317;%-0.036317
for i=2:nmax,
gt(i)=B3.*gt(i-1)+r(i);
end
 
gt=gt/max(abs(gt));
gt=gt';
%gammar0=gammar.*r;
%gammar=gammar';
%testgammar0=gt;
testgammar0=testgammar.*gt;
testgammar10=testgammar+0.001.*testgammar0;
testgammar20=testgammar+0.002.*testgammar0;
testgammar30=testgammar+0.003.*testgammar0;
testgammar40=testgammar+0.004.*testgammar0;
testgammar50=testgammar+0.005.*testgammar0;
 
fid = fopen('tinerror10.dat','w');
%testerror10=[testden10';testporo10';testgammar10'];
testerror10=[testden10';testporo10'];
fprintf(fid,'%6.4f %6.4f\n',testerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror10);
fclose(fid);
 
fid = fopen('tinerror20.dat','w');
%testerror20=[testden20';testporo20';testgammar20'];
testerror20=[testden20';testporo20'];
fprintf(fid,'%6.4f %6.4f\n',testerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror20);
fclose(fid);
 
fid = fopen('tinerror30.dat','w');
%testerror30=[testden30';testporo30';testgammar30'];
testerror30=[testden30';testporo30'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror30);
fprintf(fid,'%6.4f %6.4f\n',testerror30);
fclose(fid);
 
fid = fopen('tinerror40.dat','w');
%testerror40=[testden40';testporo40';testgammar40'];
testerror40=[testden40';testporo40'];
fprintf(fid,'%6.4f %6.4f\n',testerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror40);
fclose(fid);
 
fid = fopen('tinerror50.dat','w');
%testerror50=[testden50';testporo50';testgammar50'];
testerror50=[testden50';testporo50'];
fprintf(fid,'%6.4f %6.4f\n',testerror50);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror50);
fclose(fid);
 
load tinerror10.dat;
load tinerror20.dat;
load tinerror30.dat;
load tinerror40.dat;
load tinerror50.dat
 
 
[e10t]= tramnmx(tinerror10',minp,maxp);
[p2e10t] = trapca(e10t,transMat);
a10t = svmval(p2e10t',xsup,w,w0,kernel,kerneloption);
%a10t = gpfwd(net,p2e10t');
[a10t] = postmnmx(a10t',mint,maxt);
 
 
%tat=tat';
tdev10=tat'-a10t;
tdev10=tdev10';
tdp10=tdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
tstmse10=mse(tdp10);
 
figure(45);
 
subplot(1,1,1),plot(tdp10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e20t]= tramnmx(tinerror20',minp,maxp);
[p2e20t] = trapca(e20t,transMat);
a20t = svmval(p2e20t',xsup,w,w0,kernel,kerneloption);
%a20t = gpfwd(net,p2e20t');
[a20t] = postmnmx(a20t',mint,maxt);
 
 
 
tdev20=tat'-a20t;
tdev20=tdev20';
tdp20=tdev20(:,1);
tstmse20=mse(tdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(46);
 
subplot(1,1,1),plot(tdp20,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
[e30t]= tramnmx(tinerror30',minp,maxp);
[p2e30t] = trapca(e30t,transMat);
a30t = svmval(p2e30t',xsup,w,w0,kernel,kerneloption);
%a30t = gpfwd(net,p2e10t');
[a30t] = postmnmx(a30t',mint,maxt);
 
 
 
tdev30=tat'-a30t;
tdev30=tdev30';
tdp30=tdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
tstmse30=mse(tdp30);
 
figure(47);
 
subplot(1,1,1),plot(tdp30,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
[e40t]= tramnmx(tinerror40',minp,maxp);
[p2e40t] = trapca(e40t,transMat);
a40t = svmval(p2e40t',xsup,w,w0,kernel,kerneloption);
%a40t = gpfwd(net,p2e40t');
[a40t] = postmnmx(a40t',mint,maxt);
 
tdev40=tat'-a40t;
tdev40=tdev40';
tdp40=tdev40(:,1);
tstmse40=mse(tdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(48);
 
subplot(1,1,1),plot(tdp40,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
[e50t]= tramnmx(tinerror50',minp,maxp);
[p2e50t] = trapca(e50t,transMat);
a50t = svmval(p2e50t',xsup,w,w0,kernel,kerneloption);
%a50t = gpfwd(net,p2e50t');
[a50t] = postmnmx(a50t',mint,maxt);
 
 
 
tdev50=tat'-a50t;
tdev50=tdev50';
tdp50=tdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
tstmse50=mse(tdp50);
 
figure(50);
 
subplot(1,1,1),plot(tdp50,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
 
figure(51)
nx=[10 20 30 40 50];
msetst=[tstmse10 tstmse20 tstmse30 tstmse40 tstmse50];
bar(nx,msetst);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Test Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Obs.','SVMreg', 'SVMreg+STD','SVMreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
figure(52)
count=[msetr' mseval' msetst']
ybar=count(1:5,:);
bar(ybar,'LineWidth',3,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('MSE Analysis in Noisy Input Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
meancount=mean(count')
 
figure(53)
count=[mseval' msetst'];
ybar=count(1:5,:);
bar(ybar,'LineWidth',3,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('MSE Analysis in Noisy Input Data');
legend('Validation','Test','Location','northeast');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
plotregression(tr,a10tr')
xlabel('Target','FontName','Tahoma','Fontsize',20);ylabel('SVMreg','FontName','Tahoma','Fontsize',20);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Training Data+10% Noise','Location','northeast');
 
 
plotregression(tav,a10v')
xlabel('Target','FontName','Tahoma','Fontsize',20);ylabel('SVMreg','FontName','Tahoma','Fontsize',20);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
 
 
plotregression(tat,a10t')
xlabel('Target','FontName','Tahoma','Fontsize',20);ylabel('SVMreg','FontName','Tahoma','Fontsize',20);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Test Data+10% Noise','Location','northeast');
 
%correlation analysis at different target noise
vinerror10=[ff,gg];
[e10]= tramnmx(vinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
a10 = svmval(p2e10',xsup,w,w0,kernel,kerneloption);
%a10 = gpfwd(net,p2e10');
[a10] = postmnmx(a10',mint,maxt);
 
noisetarget10tr=target10(:,iitr);
noisetarget20tr=target20(:,iitr);
noisetarget30tr=target30(:,iitr);
noisetarget40tr=target40(:,iitr);
noisetarget50tr=target50(:,iitr);
 
noisetarget10tv=target10(:,iival);
noisetarget20tv=target20(:,iival);
noisetarget30tv=target30(:,iival);
noisetarget40tv=target40(:,iival);
noisetarget50tv=target50(:,iival);
 
noisetarget10tst=target10(:,iitst);
noisetarget20tst=target20(:,iitst);
noisetarget30tst=target30(:,iitst);
noisetarget40tst=target40(:,iitst);
noisetarget50tst=target50(:,iitst);
 
noisea10tr=a10(:,iitr);
noisea10tv=a10(:,iival);
noisea10tst=a10(:,iitst);
 
%%%%%%%%%% R analysis at different level of target noise/trmsecorrcoef10=mse(noisetarget10tr-noisea10tr);
%10
trcorrcoef10=corrcoef(noisetarget10tr,noisea10tr)
valcorrcoef10=corrcoef(noisetarget10tv,noisea10tv)
testcorrcoef10=corrcoef(noisetarget10tst,noisea10tst)
 
noisetrcod10=diag(trcorrcoef10,1).*diag(trcorrcoef10,1)
noisevalcod10=diag(valcorrcoef10,1).*diag(valcorrcoef10,1)
noisetestcod10=diag(testcorrcoef10,1).*diag(testcorrcoef10,1)
%%%%20
trcorrcoef20=corrcoef(noisetarget20tr,noisea10tr)
valcorrcoef20=corrcoef(noisetarget20tv,noisea10tv)
testcorrcoef20=corrcoef(noisetarget20tst,noisea10tst)
%30
noisetrcod20=diag(trcorrcoef20,1).*diag(trcorrcoef20,1)
noisevalcod20=diag(valcorrcoef20,1).*diag(valcorrcoef20,1)
noisetestcod20=diag(testcorrcoef20,1).*diag(testcorrcoef20,1)
%%%%%%%%%%30
trcorrcoef30=corrcoef(noisetarget30tr,noisea10tr)
valcorrcoef30=corrcoef(noisetarget30tv,noisea10tv)
testcorrcoef30=corrcoef(noisetarget30tst,noisea10tst)
%50
noisetrcod30=diag(trcorrcoef30,1).*diag(trcorrcoef30,1)
noisevalcod30=diag(valcorrcoef30,1).*diag(valcorrcoef30,1)
noisetestcod30=diag(testcorrcoef30,1).*diag(testcorrcoef30,1)
%%%%%%%%%%%%%%%%%%%40
trcorrcoef40=corrcoef(noisetarget40tr,noisea10tr)
valcorrcoef40=corrcoef(noisetarget40tv,noisea10tv)
testcorrcoef40=corrcoef(noisetarget40tst,noisea10tst)
 
noisetrcod40=diag(trcorrcoef40,1).*diag(trcorrcoef40,1)
noisevalcod40=diag(valcorrcoef40,1).*diag(valcorrcoef40,1)
noisetestcod40=diag(testcorrcoef40,1).*diag(testcorrcoef40,1)
%%%%%%50
trcorrcoef50=corrcoef(noisetarget50tr,noisea10tr)
valcorrcoef50=corrcoef(noisetarget50tv,noisea10tv)
testcorrcoef50=corrcoef(noisetarget50tst,noisea10tst)
 
noisetrcod50=diag(trcorrcoef50,1).*diag(trcorrcoef50,1)
noisevalcod50=diag(valcorrcoef50,1).*diag(valcorrcoef50,1);
noisetestcod50=diag(testcorrcoef50,1).*diag(testcorrcoef50,1)
 
r2tr=[noisetrcod10 noisetrcod20 noisetrcod30 noisetrcod40 noisetrcod50];
r2tv=[noisevalcod10 noisevalcod20 noisevalcod30 noisevalcod40 noisevalcod50];
r2tst=[noisetestcod10 noisetestcod20 noisetestcod30 noisetestcod40 noisetestcod50];
 
r2tr=sqrt(r2tr);
r2tv=sqrt(r2tv);
r2tst=sqrt(r2tst);
 
figure(70)
countr2=[r2tr' r2tv' r2tst']
ybar2=countr2(1:5,:);
bar(ybar2,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('R ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('R Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2=mean(countr2')
 
%..................MSE analysis
trmsecorrcoef10=mse(noisetarget10tr-noisea10tr);
valmsecorrcoef10=mse(noisetarget10tv-noisea10tv);
testmsecorrcoef10=mse(noisetarget10tst-noisea10tst);
 
trmsecorrcoef20=mse(noisetarget20tr-noisea10tr);
valmsecorrcoef20=mse(noisetarget20tv-noisea10tv);
testmsecorrcoef20=mse(noisetarget20tst-noisea10tst);
 
trmsecorrcoef30=mse(noisetarget30tr-noisea10tr);
valmsecorrcoef30=mse(noisetarget30tv-noisea10tv);
testmsecorrcoef30=mse(noisetarget30tst-noisea10tst);
 
trmsecorrcoef40=mse(noisetarget40tr-noisea10tr);
valmsecorrcoef40=mse(noisetarget40tv-noisea10tv);
testmsecorrcoef40=mse(noisetarget40tst-noisea10tst);
 
trmsecorrcoef50=mse(noisetarget50tr-noisea10tr);
valmsecorrcoef50=mse(noisetarget50tv-noisea10tv);
testmsecorrcoef50=mse(noisetarget50tst-noisea10tst);
 
r2msetr=[trmsecorrcoef10 trmsecorrcoef20 trmsecorrcoef30 trmsecorrcoef40 trmsecorrcoef50];
r2msetv=[valmsecorrcoef10 valmsecorrcoef20 valmsecorrcoef30 valmsecorrcoef40 valmsecorrcoef50];
r2msetst=[testmsecorrcoef10 testmsecorrcoef20 testmsecorrcoef30 testmsecorrcoef40 testmsecorrcoef50];
 
figure(71)
countr2mse=[r2msetr' r2msetv' r2msetst'];
ybar2mse=countr2mse(1:5,:);
bar(ybar2mse,'LineWidth',4,'BarWidth',1);
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('MSE Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2mse=mean(countr2mse')
 
%Input noise analysis based on R analysis
%%%%%%%%%% R analysis
in10tr=a10tr;
in20tr=a20tr;
in30tr=a30tr;
in40tr=a40tr;
in50tr=a50tr;
 
in10tv=a10v;
in20tv=a20v;
in30tv=a30v;
in40tv=a40v;
in50tv=a50v;
 
in10tst=a10t;
in20tst=a20t;
in30tst=a30t;
in40tst=a40t;
in50tst=a50t;
 
out10tr=tr;
out20tr=tr;
out30tr=tr;
out40tr=tr;
out50tr=tr;
 
out10tv=tav;
out20tv=tav;
out30tv=tav;
out40tv=tav;
out50tv=tav;
 
out10tst=tat;
out20tst=tat;
out30tst=tat;
out40tst=tat;
out50tst=tat;
 
trcorrcoef10=corrcoef(in10tr,out10tr);
valcorrcoef10=corrcoef(in10tv,out10tv);
testcorrcoef10=corrcoef(in10tst,out10tst);
 
noisetrcod10=diag(trcorrcoef10,1).*diag(trcorrcoef10,1);
noisevalcod10=diag(valcorrcoef10,1).*diag(valcorrcoef10,1);
noisetestcod10=diag(testcorrcoef10,1).*diag(testcorrcoef10,1);
%%%%
trcorrcoef20=corrcoef(in20tr,out20tr);
valcorrcoef20=corrcoef(in20tv,out20tv);
testcorrcoef20=corrcoef(in20tst,out20tst);
 
noisetrcod20=diag(trcorrcoef20,1).*diag(trcorrcoef20,1);
noisevalcod20=diag(valcorrcoef20,1).*diag(valcorrcoef20,1);
noisetestcod20=diag(testcorrcoef20,1).*diag(testcorrcoef20,1);
%%%%%%%%%%
trcorrcoef30=corrcoef(in30tr,out30tr);
valcorrcoef30=corrcoef(in30tv,out30tv);
testcorrcoef30=corrcoef(in30tst,out30tst);
 
noisetrcod30=diag(trcorrcoef30,1).*diag(trcorrcoef30,1);
noisevalcod30=diag(valcorrcoef30,1).*diag(valcorrcoef30,1);
noisetestcod30=diag(testcorrcoef30,1).*diag(testcorrcoef30,1);
%%%%%%%%%%%%%%%%%%%5
trcorrcoef40=corrcoef(in40tr,out40tr);
valcorrcoef40=corrcoef(in40tv,out40tv);
testcorrcoef40=corrcoef(in40tst,out40tst);
 
noisetrcod40=diag(trcorrcoef40,1).*diag(trcorrcoef40,1);
noisevalcod40=diag(valcorrcoef40,1).*diag(valcorrcoef40,1);
noisetestcod40=diag(testcorrcoef40,1).*diag(testcorrcoef40,1);
%%%%%%
trcorrcoef50=corrcoef(in50tr,out50tr);
valcorrcoef50=corrcoef(in50tv,out50tv);
testcorrcoef50=corrcoef(in50tst,out50tst);
 
noisetrcod50=diag(trcorrcoef50,1).*diag(trcorrcoef50,1);
noisevalcod50=diag(valcorrcoef50,1).*diag(valcorrcoef50,1);
noisetestcod50=diag(testcorrcoef50,1).*diag(testcorrcoef50,1);
 
r2tr=[noisetrcod10 noisetrcod20 noisetrcod30 noisetrcod40 noisetrcod50];
r2tv=[noisevalcod10 noisevalcod20 noisevalcod30 noisevalcod40 noisevalcod50];
r2tst=[noisetestcod10 noisetestcod20 noisetestcod30 noisetestcod40 noisetestcod50];
 
r2tr=sqrt(r2tr);
r2tv=sqrt(r2tv);
r2tst=sqrt(r2tst);
 
figure(72)
countr2=[r2tr' r2tv' r2tst'];
ybar2=countr2(1:5,:);
bar(ybar2,'LineWidth',4,'BarWidth',1);
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('R ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('R Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2=mean(countr2')
 
%equation for calculating alphasm
bias=0.3346;
lx=aa(:,iitr);
ly=bb(:,iitr);
sigsm=std(tr);
twosigmasq=2.*sigsm.*sigsm;
distance=(lx-mean(lx))*(ly-mean(ly))';
distance=distance./twosigmasq;
alphasm=(tr-bias)./exp(-distance);
 
figure(73)
plot(alphasm,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
%trp0=alphasm.*exp(-distance)+bias;
% plot(trp0);hold on
% plot(tr)
xlabel('Training dataset','FontName','Tahoma','Fontsize',24);ylabel('alpha ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
%title('R Analysis in Noisy Target Data');
%legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
%Training data;
totalinerror=totalinerror10;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
aan=aan';
bbn=bbn';
 
bias=0.3098;
lx=aan(:,iitr);
ly=bbn(:,iitr);
sigsm10=std(target10);
twosigmasq10=2.*sigsm10.*sigsm10;
distance=(lx-mean(lx)).*(ly-mean(ly));
distance=distance./twosigmasq10;
alphasm10=(target10(:,iitr)-bias)./exp(-distance);
 
%for 20%
totalinerror=totalinerror20;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
aan=aan';
bbn=bbn';
bias=0.3098;
lx=aan(:,iitr);
ly=bbn(:,iitr);
sigsm20=std(target20);
twosigmasq20=2.*sigsm20.*sigsm20;
distance=(lx-mean(lx)).*(ly-mean(ly));
distance=distance./twosigmasq20;
alphasm20=(target20(:,iitr)-bias)./exp(-distance);
%for 30%
totalinerror=totalinerror30;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
aan=aan';
bbn=bbn';
bias=0.3098;
lx=aan(:,iitr);
ly=bbn(:,iitr);
sigsm30=std(target30);
twosigmasq30=2.*sigsm30.*sigsm30;
distance=(lx-mean(lx)).*(ly-mean(ly));
distance=distance./twosigmasq30;
alphasm30=(target30(:,iitr)-bias)./exp(-distance);
 
%for 40%
totalinerror=totalinerror40;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
aan=aan';
bbn=bbn';
bias=0.3098;
lx=aan(:,iitr);
ly=bbn(:,iitr);
sigsm40=std(target40);
twosigmasq40=2.*sigsm40.*sigsm40;
distance=(lx-mean(lx)).*(ly-mean(ly));
distance=distance./twosigmasq40;
alphasm40=(target40(:,iitr)-bias)./exp(-distance);
%for 50%
totalinerror=totalinerror50;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
aan=aan';
bbn=bbn';
bias=0.3098;
lx=aan(:,iitr);
ly=bbn(:,iitr);
sigsm50=std(target50);
twosigmasq50=2.*sigsm50.*sigsm50;
distance=(lx-mean(lx)).*(ly-mean(ly));
distance=distance./twosigmasq50;
alphasm50=(target50(:,iitr)-bias)./exp(-distance);
 
figure(74) 
plot(alphasm10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
            plot(alphasm20,'bo-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
            plot(alphasm30,'go-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
            plot(alphasm40,'co-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
            plot(alphasm50,'mo-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
%trp50=alphasm50.*exp(-distance)+bias;
% plot(trp50);hold on
% plot(tr)
 
xlabel('Training dataset','FontName','Tahoma','Fontsize',24);ylabel('alpha ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
%title('R Analysis in Noisy Target Data');
%legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight

 
% %DEMEV1   Demonstrate Bayesian regression for the MLP...............
% % Gaussian Process 
% % Description
% 
% 
clc;
clear all;
close all;
 
%Training data;
load gravdeep.dat;
aa=gravdeep(:,1);%longitude
bb=gravdeep(:,2);%latitude
cc=gravdeep(:,3);%BG value
 
dd=gravdeep(:,4);%Altitude
ee=gravdeep(:,5);%deep basement
% %----------------------------------------------------
 
totalden = aa;
totalporo = bb;
totalgammar = cc;
%........................................
 
%Noise stability test
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
Aa1=0.99;%0.99
dsm(1)=a(1);
for i=2:nmax,
dsm(i)=Aa1.*dsm(i-1)+r(i);
end
dsm=dsm/max(abs(dsm));
dsm=dsm';
%r=r';
%den0=den.*r;
%den=a;
%totalden0=dsm;
totalden0=totalden.*dsm;
totalden10=totalden+0.001.*totalden0;
totalden20=totalden+0.002.*totalden0;
totalden30=totalden+0.003.*totalden0;
totalden40=totalden+0.004.*totalden0;
totalden50=totalden+0.005.*totalden0;
 
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
pr(1)=r(1);
Ba2=-0.0957;%-0.0957
for i=2:nmax,
pr(i)=Ba2.*pr(i-1)+r(i);
end
 
pr=pr/max(abs(pr));
pr=pr';
%poro0=poro.*r;
%poro=b;
%totalporo0=pr;
totalporo0=totalporo.*pr;
totalporo10=totalporo+0.001.*totalporo0;
totalporo20=totalporo+0.002.*totalporo0;
totalporo30=totalporo+0.003.*totalporo0;
totalporo40=totalporo+0.004.*totalporo0;
totalporo50=totalporo+0.005.*totalporo0;
 
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
gg(1)=r(1);
Ca3=0.0142;%0.0142
for i=2:nmax,
gg(i)=Ca3.*gg(i-1)+r(i);
end
gg=gg/max(abs(gg));
gg=gg';
%gammar0=gammar.*r;
%gammar=gammar';
%totalgammar0=gg;
totalgammar0=totalgammar.*gg;
totalgammar10=totalgammar+0.001.*totalgammar0;
totalgammar20=totalgammar+0.002.*totalgammar0;
totalgammar30=totalgammar+0.003.*totalgammar0;
totalgammar40=totalgammar+0.004.*totalgammar0;
totalgammar50=totalgammar+0.005.*totalgammar0;
 
fid = fopen('totalinerror10.dat','w');
totalerror10=[totalden10';totalporo10';totalgammar10'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror10);
fclose(fid);
 
fid = fopen('totalinerror20.dat','w');
totalerror20=[totalden20';totalporo20';totalgammar20'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror20);
fclose(fid);
 
fid = fopen('totalinerror30.dat','w');
totalerror30=[totalden30';totalporo30';totalgammar30'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror30);
fclose(fid);
 
fid = fopen('totalinerror40.dat','w');
totalerror40=[totalden40';totalporo40';totalgammar40'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror40);
fclose(fid);
 
fid = fopen('totalinerror50.dat','w');
totalerror50=[totalden50';totalporo50';totalgammar50'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror50);
fclose(fid);
 
load totalinerror10.dat
load totalinerror20.dat
load totalinerror30.dat
load totalinerror40.dat
load totalinerror50.dat
 
%Training data;
totalinerror=totalinerror30;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
ccn=totalinerror(:,3);%BG value
 
%........................................................................
% Real data
load gravwhole.dat;
ff=gravwhole(:,1);%longitude
gg=gravwhole(:,2);%latitude
hh=gravwhole(:,3);%BG value
ii=gravwhole(:,4);%altitude
 
x=[aan,bbn];
t=ee;
 
x=x';
t=t';
 
%.....................................................noise adding......
 %ndata=211;
 noise = 0.15;          % Standard deviation of noise distribution.
 
 %*********************************************************
target=t;
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
Ta1=0.99;%0.99
dm(1)=a(1);
for i=2:nmax,
dm(i)=Ta1.*dm(i-1)+r(i);
end
 
dm=dm/max(abs(dm));
target0=target.*dm;
target10=target+0.1.*target0;
target20=target+0.2.*target0;
target30=target+0.3.*target0;
target40=target+0.4.*target0;
target50=target+0.5.*target0;
t=target30;
 
%**********************************************************************
 
%***********************************************************************
[xn,minp,maxp,tn,mint,maxt] = premnmx(x,t);
[xtrans,transMat] = prepca(xn,0.0000000000000002);
[R,Q] = size(xtrans);
iitst = 2:4:Q;
iival = 4:4:Q;
iitr = [1:4:Q 3:4:Q];
valP = xtrans(:,iival); val.T = tn(:,iival);
testP = xtrans(:,iitst); test.T = tn(:,iitst);
xtr = xtrans(:,iitr); ttr = tn(:,iitr);
xr=xtr';
zr=ttr';
 
 
 
net = gp(2, 'sqexp');
%net = gp(2, 'ratquad');
% Initialise the parameters.
prior.pr_mean = 0;
prior.pr_var = 0.01;
net = gpinit(net, xr, zr, prior);
 
% Now train to find the hyperparameters.
options = foptions;
options(1) = 1;
options(14) = 300;
 
[net, options] = netopt(net, options, xr, zr, 'scg');
 
rel = exp(net.inweights);
 
 
[an, sqxr] = gpfwd(net,xr);
[y] = postmnmx(an',mint,maxt);
sqtr=2.00*sqrt(sqxr');
ee=ee';
tr=ee(:,iitr);
 
 
figure(26);
 
xis=1:106;
E1=sqtr(:,1).*ones(size(y));
%subplot(3,1,1),
errorbar(xis,y,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tr,'ro-',xis, y,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -6 6]);
legend('Error-bar','Target','GPreg.');
xlabel('No. of Data','Fontsize',20);
ylabel('Sediment Depth(Km)','Fontsize',20);
title('Training Interval')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%........................................
figure;
an = gpfwd(net,xr);
[y] = postmnmx(an',mint,maxt);
msigtr=sqtr;
%ee=ee';
tr=t(:,iitr);
plot(tr,'-r','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
plot(y,'--b>','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
 plot(y + msigtr, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
 plot(y - msigtr, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
xlabel('Number of Training Samples','FontName','Tahoma','Fontsize',20);ylabel('Shallow Interface Depth(Km)','FontName','Tahoma','Fontsize',20);
title('Shallow Interface Depth Estimation:Training Interval','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
legend('Obs.','GPreg', 'GPreg+STD','GPreg-STD');
 
plotregression(tr,y)
xlabel('Obs.','FontName','Tahoma','Fontsize',20);ylabel('GPreg','FontName','Tahoma','Fontsize',20);
%title('Basement Depth Prediction:Training Interval','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);
 
%..........................................

 
%............................................
%sqvalx=sigval;
[av,sqvalx] = gpfwd(net,valP');
[avy] = postmnmx(av',mint,maxt);
sqval=2.00*sqrt(sqvalx');
[av] = postmnmx(av',mint,maxt);
%ee=ee';
tav=ee(:,iival);
 
figure(27);
xis=1:52;
E2=sqtr(:,1).*ones(size(av));
 
%subplot(3,1,1),
errorbar(xis,avy,E2,'c');grid on;axis([0 52 0 6]);hold on
%title('Errobar plot');hold on
plot(xis,tav,'ro-',xis, avy,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 52 0 6]);
legend('Error-bar','Target','GPreg.');
%plot(tat);hold on;plot(ate);
xlabel('No. of Data','Fontsize',20);
ylabel('Sediment Depth(Km)','Fontsize',20);
title('Validation Interval')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
msigval=sqval;
 
figure;
av = gpfwd(net,valP');
[av] = postmnmx(av',mint,maxt);
%tav=ee(:,iival);
plot(tav,'-r','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
plot(av,'--b>','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
 plot(av + msigval, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
 plot(av - msigval, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
xlabel('Number of Validation Samples','FontName','Tahoma','Fontsize',20);ylabel('Shallow Interface Depth(Km) ','FontName','Tahoma','Fontsize',20);
title('Shallow Interface Depth Estimation:Validation Interval','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
legend('Obs.','GPreg', 'GPreg+STD','GPreg-STD');
 
 
plotregression(tav,av)
xlabel('Obs.','FontName','Tahoma','Fontsize',20);ylabel('GPreg','FontName','Tahoma','Fontsize',20);
%title('Basement Depth Prediction:Training Interval','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);
 
%sqtestx=sigtst;
[ate,sqtestx] = gpfwd(net,testP');
[atey] = postmnmx(ate',mint,maxt);
sqtest=2.00*sqrt(sqtestx);
msigtst=sqtest;
%tat=t(:,iitst);
tat=ee(:,iitst);
 
figure(28);
xis=1:53;
E3=sqtr(:,1).*ones(size(atey));
 
%subplot(3,1,1),
errorbar(xis,atey,E3,'c');grid on;axis([0 53 0 6]);hold on
%title('Errobar plot');hold on
plot(xis,tat,'ro-',xis, atey,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 53 0 6]);
legend('Error-bar','Target','GPreg.');
%plot(tat);hold on;plot(ate);
xlabel('No. of Data','Fontsize',20);
ylabel('Sediment Depth(Km)','Fontsize',20);
title('Test Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
figure
ate = gpfwd(net,testP');
[ate] = postmnmx(ate',mint,maxt);
%tat=ee(:,iitst);
plot(tat,'-r','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
plot(ate,'--b>','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6); hold on;
 plot(ate + msigtst, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
 plot(ate - msigtst, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
xlabel('Number of test Samples','FontName','Tahoma','Fontsize',20);ylabel('Shallow Interface Depth(Km) ','FontName','Tahoma','Fontsize',20);
title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);
legend('Obs.','GPreg', 'GPreg+STD','GPreg-STD');
 
 
plotregression(tat,ate)
xlabel('Obs.','FontName','Tahoma','Fontsize',20);ylabel('GPreg','FontName','Tahoma','Fontsize',20);
%title('Basement Depth Prediction:Training Interval','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);
 
figure
% Real data analysis;
vinerror10=[ff,gg];
[e10]= tramnmx(vinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
a10 = gpfwd(net,p2e10');
[a10] = postmnmx(a10',mint,maxt);
%[a10]=-a10;
 
fdata=488;
gdata=488;
fmax=max(ff);
gmax=max(gg);
fmin=min(ff);
gmin=min(gg);
fv = linspace(fmin, fmax,fdata);
gv = linspace(gmin, gmax,gdata);
[fi,gi] = meshgrid(fv,gv);
c1i = griddata(ff,gg,a10,fi,gi,'cubic');
 
meshz(fi,gi,c1i);hold on;
colorbar;
plot3(ff,gg,a10,'v','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
%view(10,90);
% view(0,90);
%shading interp;
grid on;
box on;
%axis tight;
%axis auto;
%axis square
colorbar;
set(gca,'LineWidth',2,'FontName','Tahoma','Fontsize',20);axis tight;
xlabel('Longitude(in degree)','Fontsize',20);
ylabel('Latitude(in degree)','Fontsize',20);
zlabel('Deeper Interface Depth(Km)','Fontsize',20);set(gca,'Zdir','reverse');
title('Baement Depth Estimation','Fontsize',20);
%.........................................................................
load datav2.dat;
sigmain1=datav2(:,1);
 
c2i = griddata(ff,gg,sigmain1,fi,gi,'cubic');
 
meshz(fi,gi,c2i);hold on;%meshc(fi,gi,c1i);hold on;
plot3(ff,gg,sigmain1,'v','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
%view(10,90);
% view(0,90);
%shading interp;
grid on;
box on;
%axis tight;
%axis auto;
%axis square
colorbar;
colormap(jet);
set(gca,'LineWidth',2,'FontName','Tahoma','Fontsize',20);axis tight;
xlabel('Longitude(in degree)','Fontsize',20);
ylabel('Latitude(in degree)','Fontsize',20);
zlabel('STD (Km)','Fontsize',20);%set(gca,'Zdir','reverse');
title('Uncertainty Estimation in Deeper Interface Depth Prediction ','Fontsize',20);
%vt = t(:,iival);
 
 
%clear all
%Statistical %Training....................................................
%.........................................................................
ntotal=211;
ntrain=106;
nval=52;
ntest=53;
 
x1train=tr(1:ntrain);
x1val=tav(1:nval);
x1test=tat(1:ntest);
 
x2train=y(1:ntrain);
x2val=av(1:nval);
x2test=ate(1:ntest);
 
disp('Statistics of x1/obs......................................');
x1avgtrain=mean(x1train)
x1avgnval=mean(x1val)
x1avgtest=mean(x1test)
 
x1stdtrain=std(x1train)
x1stdval=std(x1val)
x1stdtest=std(x1test)
 
x1trmedian=median(x1train)
x1valmedian=median(x1val)
x1testmedian=median(x1test)
 
x1trmode=mode(x1train)
x1valmode=mode(x1val)
x1testmode=mode(x1test)
 
x1trskewness=skewness(x1train)
x1valskewness=skewness(x1val)
x1testskewness=skewness(x1test)
 
x1trkurtosis=kurtosis(x1train)
x1valkurtosis=kurtosis(x1val)
x1testkurtosis=kurtosis(x1test)
 
disp('Statistics of x2/GPreg......................................')
x2avgtrain=mean(x2train)
x2avgnval=mean(x2val)
x2avgtest=mean(x2test)
 
x2stdtrain=std(x2train)
x2stdval=std(x2val)
x2stdtest=std(x2test)
 
x2trmedian=median(x2train);
x2valmedian=median(x2val)
x2testmedian=median(x2test)
 
x2trmode=mode(x2train)
x2valmode=mode(x2val)
x2testmode=mode(x2test)
 
x2trskewness=skewness(x2train)
x2valskewness=skewness(x2val)
x2testskewness=skewness(x2test)
 
x2trkurtosis=kurtosis(x2train)
x2valkurtosis=kurtosis(x2val)
x2testkurtosis=kurtosis(x2test)
 
x2errtrain=mse(x2train-x1train)
x2errval=mse(x2val-x1val)
x2errtest=mse(x2test-x1test)
 
x2maetrain=mae(x2train-x1train)
x2maeval=mae(x2val-x1val)
x2maetest=mae(x2test-x1test)
 
x2retrain=1.0-sum((x1train-x2train).^2)./sum(x1train.^2)
x2reval=1.0-sum((x1val-x2val).^2)./sum(x1val.^2);
x2retest=1.0-sum((x1test-x2test).^2)./sum(x1test.^2)
 
x2dtrain=1.0-sum((x1train-x2train).^2)./sum(abs(x2train-mean(x1train))+abs(x1train-mean(x1train)).^2)
x2dval=1.0-sum((x1val-x2val).^2)./sum(abs(x2val-mean(x1val))+abs(x1val-mean(x1val)).^2)
x2dtest=1.0-sum((x1test-x2test).^2)./sum(abs(x2test-mean(x1test))+abs(x1test-mean(x1test)).^2)
 
trcorrcoef=corrcoef(x1train,x2train)
valcorrcoef=corrcoef(x1val,x2val)
testcorrcoef=corrcoef(x1test,x2test)
 
trcod=diag(trcorrcoef,1).*diag(trcorrcoef,1)
valcod=diag(valcorrcoef,1).*diag(valcorrcoef,1)
testcod=diag(testcorrcoef,1).*diag(testcorrcoef,1)
 
%***************Training Error Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error analysis**********************************************************************
aa=aa';
bb=bb';
cc=cc';
 
inputden = aa(:,iitr);
inputporo = bb(:,iitr);
inputgammar = cc(:,iitr);
 
inputden=inputden';
inputporo=inputporo';
inputgammar=inputgammar';
 
randn('state',0);
nmax=106;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
a1=0.94;%0.94
dtr(1)=a(1);
for i=2:nmax,
dtr(i)=a1.*dtr(i-1)+r(i);
end
 
dtr=dtr/max(abs(dtr));
dtr=dtr';
%r=r';
%den0=den.*r;
%den=a;
%inputden0=dtr;
inputden0=inputden.*dtr;
inputden10=inputden+0.001.*inputden0;
inputden20=inputden+0.002.*inputden0;
inputden30=inputden+0.003.*inputden0;
inputden40=inputden+0.004.*inputden0;
inputden50=inputden+0.005.*inputden0;
 
randn('state',0);
nmax=106;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
ptr(1)=r(1);
a2=-0.107;%-0.107
for i=2:nmax,
ptr(i)=a2.*ptr(i-1)+r(i);
end
ptr=ptr';
ptr=ptr/max(abs(ptr));
%poro0=poro.*r;
%poro=b;
%inputporo0=ptr;
inputporo0=inputporo.*ptr;
inputporo10=inputporo+0.001.*inputporo0;
inputporo20=inputporo+0.002.*inputporo0;
inputporo30=inputporo+0.003.*inputporo0;
inputporo40=inputporo+0.004.*inputporo0;
inputporo50=inputporo+0.005.*inputporo0;
 
randn('state',0);
nmax=106;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
ggtr(1)=r(1);
a3=0.122;%0.122
for i=2:nmax,
ggtr(i)=a3.*ggtr(i-1)+r(i);
end
 
ggtr=ggtr/max(abs(ggtr));
ggtr=ggtr';
%gammar0=gammar.*r;
%gammar=gammar';
%inputgammar0=ggtr;
inputgammar0=inputgammar.*ggtr;
inputgammar10=inputgammar+0.001.*inputgammar0;
inputgammar20=inputgammar+0.002.*inputgammar0;
inputgammar30=inputgammar+0.003.*inputgammar0;
inputgammar40=inputgammar+0.004.*inputgammar0;
inputgammar50=inputgammar+0.005.*inputgammar0;
 
fid = fopen('sinerror10.dat','w');
%inputerror10=[inputden10';inputporo10';inputgammar10'];
inputerror10=[inputden10';inputporo10'];
fprintf(fid,'%6.4f %6.4f\n',inputerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror10);
fclose(fid);
 
fid = fopen('sinerror20.dat','w');
%inputerror20=[inputden20';inputporo20';inputgammar20'];
inputerror20=[inputden20';inputporo20'];
fprintf(fid,'%6.4f %6.4f\n',inputerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror20);
fclose(fid);
 
fid = fopen('sinerror30.dat','w');
%inputerror30=[inputden30';inputporo30';inputgammar30'];
inputerror30=[inputden30';inputporo30'];
fprintf(fid,'%6.4f %6.4f\n',inputerror30);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror30);
fclose(fid);
 
fid = fopen('sinerror40.dat','w');
%inputerror40=[inputden40';inputporo40';inputgammar40'];
inputerror40=[inputden40';inputporo40'];
fprintf(fid,'%6.4f %6.4f\n',inputerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror40);
fclose(fid);
 
fid = fopen('sinerror50.dat','w');
%inputerror50=[inputden50';inputporo50';inputgammar50'];
inputerror50=[inputden50';inputporo50'];
fprintf(fid,'%6.4f %6.4f\n',inputerror50);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror50);
fclose(fid);
 
load sinerror10.dat;
load sinerror20.dat;
load sinerror30.dat;
load sinerror40.dat;
load sinerror50.dat
 
 
[e10]= tramnmx(sinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
a10 = gpfwd(net,p2e10');
[a10tr] = postmnmx(a10',mint,maxt);
 
 
 
%t=t(:,iitr);
%tr=t(:,iitr);
trdev10=tr-a10tr;
trdev10=trdev10';
trdp10=trdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
trmse10=mse(trdp10);
 
figure(29);
 
subplot(1,1,1),plot(trdp10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e20tr]= tramnmx(sinerror20',minp,maxp);
[p2e20tr] = trapca(e20tr,transMat);
a20tr = gpfwd(net,p2e20tr');
[a20tr] = postmnmx(a20tr',mint,maxt);
 
 
 
trdev20=tr-a20tr;
trdev20=trdev20';
trdp20=trdev20(:,1);
trmse20=mse(trdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(30);
 
subplot(1,1,1),plot(trdp20,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e30tr]= tramnmx(sinerror30',minp,maxp);
[p2e30tr] = trapca(e30tr,transMat);
a30tr = gpfwd(net,p2e30tr');
[a30tr] = postmnmx(a30tr',mint,maxt);
 
 
 
trdev30=tr-a30tr;
trdev30=trdev30';
trdp30=trdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
trmse30=mse(trdp30);
 
figure(31);
 
subplot(1,1,1),plot(trdp30,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e40tr]= tramnmx(sinerror40',minp,maxp);
[p2e40tr] = trapca(e40tr,transMat);
a40tr = gpfwd(net,p2e40tr');
[a40tr] = postmnmx(a40tr',mint,maxt);
 
trdev40=tr-a40tr;
trdev40=trdev40';
trdp40=trdev40(:,1);
trmse40=mse(trdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(32);
 
subplot(1,1,1),plot(trdp10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e50tr]= tramnmx(sinerror50',minp,maxp);
[p2e50tr] = trapca(e50tr,transMat);
a50tr = gpfwd(net,p2e50tr');
[a50tr] = postmnmx(a50tr',mint,maxt);
 
 
 
trdev50=tr-a50tr;
trdev50=trdev50';
trdp50=trdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
trmse50=mse(trdp50);
 
figure(33);
subplot(1,1,1),plot(trdp50,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%error diagram
figure(34)
 
nx=[10 20 30 40 50];
msetr=[trmse10 trmse20 trmse30 trmse40 trmse50];
bar(nx,msetr);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Training Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Obs.','GPreg', 'GPreg+STD','GPreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%******************Validation Error ANALYSIS****************************
 
 
vdenv = aa(:,iival);
vporov = bb(:,iival);
vgammarv = cc(:,iival);
 
valden=vdenv';
valporo=vporov';
valgammar=vgammarv';
 
randn('state',0);
nmax=52;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
dv(1)=r(1);
A1=0.99;%0.99
for i=2:nmax,
dv(i)=A1.*dv(i-1)+r(i);
end
 
dv=dv/max(abs(dv));
dv=dv';
%r=r';
%den0=den.*r;
%valden0=dv;
valden0=valden.*dv;
valden10=valden+0.001.*valden0;
valden20=valden+0.002.*valden0;
valden30=valden+0.003.*valden0;
valden40=valden+0.004.*valden0;
valden50=valden+0.005.*valden0;
 
randn('state',0);
nmax=52;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
pv(1)=r(1);
A2=-0.0495;%-0.0495
for i=2:nmax,
pv(i)=A2.*pv(i-1)+r(i);
end
pv=pv';
pv=pv/max(abs(pv));
%poro0=poro.*r;
%poro=b;
%valporo0=pv;
valporo0=valporo.*pv;
valporo10=valporo+0.001.*valporo0;
valporo20=valporo+0.002.*valporo0;
valporo30=valporo+0.003.*valporo0;
valporo40=valporo+0.004.*valporo0;
valporo50=valporo+0.005.*valporo0;
 
randn('state',0);
nmax=52;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
gva(1)=r(1);
A3=0.1408;%0.1408
for i=2:nmax,
gva(i)=A3.*gva(i-1)+r(i);
end
 
gva=gva/max(abs(gva));
gva=gva';
%gammar0=gammar.*r;
%gammar=gammar';
%valgammar0=gva;
valgammar0=valgammar.*gva;
valgammar10=valgammar+0.001.*valgammar0;
valgammar20=valgammar+0.002.*valgammar0;
valgammar30=valgammar+0.003.*valgammar0;
valgammar40=valgammar+0.004.*valgammar0;
valgammar50=valgammar+0.005.*valgammar0;
 
fid = fopen('vinerror10.dat','w');
%valerror10=[valden10';valporo10';valgammar10'];
valerror10=[valden10';valporo10'];
fprintf(fid,'%6.4f %6.4f\n',valerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror10);
fclose(fid);
 
fid = fopen('vinerror20.dat','w');
%valerror20=[valden20';valporo20';valgammar20'];
valerror20=[valden20';valporo20'];
fprintf(fid,'%6.4f %6.4f\n',valerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror20);
fclose(fid);
 
fid = fopen('vinerror30.dat','w');
%valerror30=[valden30';valporo30';valgammar30'];
valerror30=[valden30';valporo30'];
fprintf(fid,'%6.4f %6.4f\n',valerror30);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror30);
fclose(fid);
 
fid = fopen('vinerror40.dat','w');
%valerror40=[valden40';valporo40';valgammar40'];
valerror40=[valden40';valporo40'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror40);
fprintf(fid,'%6.4f %6.4f\n',valerror40);
fclose(fid);
 
fid = fopen('vinerror50.dat','w');
%valerror50=[valden50';valporo50';valgammar50'];
valerror50=[valden50';valporo50'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror50);
fprintf(fid,'%6.4f %6.4f\n',valerror50);
fclose(fid);
 
load vinerror10.dat;
load vinerror20.dat;
load vinerror30.dat;
load vinerror40.dat;
load vinerror50.dat
 
 
[e10v]= tramnmx(vinerror10',minp,maxp);
[p2e10v] = trapca(e10v,transMat);
a10v = gpfwd(net,p2e10v');
[a10v] = postmnmx(a10v',mint,maxt);
 
tav=tav';
vdev10=tav'-a10v;
vdev10=vdev10';
vdp10=vdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
valmse10=mse(vdp10);
 
figure(35);
 
subplot(1,1,1),plot(vdp10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e20v]= tramnmx(vinerror20',minp,maxp);
[p2e20v] = trapca(e20v,transMat);
a20v = gpfwd(net,p2e20v');
[a20v] = postmnmx(a20v',mint,maxt);
 
vdev20=tav'-a20v;
vdev20=vdev20';
vdp20=vdev20(:,1);
valmse20=mse(vdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(36);
 
subplot(1,1,1),plot(vdp20,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e30v]= tramnmx(vinerror30',minp,maxp);
[p2e30v] = trapca(e30v,transMat);
a30v = gpfwd(net,p2e30v');
[a30v] = postmnmx(a30v',mint,maxt);
 
 
 
vdev30=tav'-a30v;
vdev30=vdev30';
vdp30=vdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
valmse30=mse(vdp30);
 
figure(37);
 
subplot(1,1,1),plot(vdp30,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e40v]= tramnmx(vinerror40',minp,maxp);
[p2e40v] = trapca(e40v,transMat);
a40v = gpfwd(net,p2e40v');
[a40v] = postmnmx(a40v',mint,maxt);
 
vdev40=tav'-a40v;
vdev40=vdev40';
vdp40=vdev40(:,1);
valmse40=mse(vdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(38);
 
subplot(1,1,1),plot(vdp40,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e50v]= tramnmx(vinerror50',minp,maxp);
[p2e50v] = trapca(e50v,transMat);
a50v = gpfwd(net,p2e50v');
[a50v] = postmnmx(a50v',mint,maxt);
 
 
 
vdev50=tav'-a50v;
vdev50=vdev50';
vdp50=vdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
valmse50=mse(vdp50);
 
figure(40);
subplot(1,1,1),plot(vdp50,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%error diagram
figure(41)
nx=[10 20 30 40 50];
mseval=[valmse10 valmse20 valmse30 valmse40 valmse50];
bar(nx,mseval);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Validation Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Obs.','GPreg', 'GPreg+STD','GPreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%**********************Test Error Analysis*****************************
%Error Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error analysis
dent = aa(:,iitst);
porot = bb(:,iitst);
gammart = cc(:,iitst);
 
testden=dent';
testporo=porot';
testgammar=gammart';
 
randn('state',0);
nmax=53;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
dtst(1)=r(1);
B1=0.99;%0.99
 
for i=2:nmax,
dtst(i)=B1.*dtst(i-1)+r(i);
end
 
dtst=dtst/max(abs(dtst));
dtst=dtst';
%r=r';
%den0=den.*r;
%den=a;
%testden0=dtst;
testden0=testden.*dtst;
testden10=testden+0.001.*testden0;
testden20=testden+0.002.*testden0;
testden30=testden+0.003.*testden0;
testden40=testden+0.004.*testden0;
testden50=testden+0.005.*testden0;
 
randn('state',0);
nmax=53;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
pt(1)=r(1);
B2=-0.2317;%-0.2317
for i=2:nmax,
pt(i)=B2.*pt(i-1)+r(i);
end
 
pt=pt/max(abs(pt));
pt=pt';
%poro0=poro.*r;
%poro=b;
%testporo0=pt;
testporo0=testporo.*pt;
testporo10=testporo+0.001.*testporo0;
testporo20=testporo+0.002.*testporo0;
testporo30=testporo+0.003.*testporo0;
testporo40=testporo+0.004.*testporo0;
testporo50=testporo+0.005.*testporo0;
 
randn('state',0);
nmax=53;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
gt(1)=r(1);
B3=-0.036317;%-0.036317
for i=2:nmax,
gt(i)=B3.*gt(i-1)+r(i);
end
 
gt=gt/max(abs(gt));
gt=gt';
%gammar0=gammar.*r;
%gammar=gammar';
%testgammar0=gt;
testgammar0=testgammar.*gt;
testgammar10=testgammar+0.001.*testgammar0;
testgammar20=testgammar+0.002.*testgammar0;
testgammar30=testgammar+0.003.*testgammar0;
testgammar40=testgammar+0.004.*testgammar0;
testgammar50=testgammar+0.005.*testgammar0;
 
fid = fopen('tinerror10.dat','w');
%testerror10=[testden10';testporo10';testgammar10'];
testerror10=[testden10';testporo10'];
fprintf(fid,'%6.4f %6.4f\n',testerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror10);
fclose(fid);
 
fid = fopen('tinerror20.dat','w');
%testerror20=[testden20';testporo20';testgammar20'];
testerror20=[testden20';testporo20'];
fprintf(fid,'%6.4f %6.4f\n',testerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror20);
fclose(fid);
 
fid = fopen('tinerror30.dat','w');
%testerror30=[testden30';testporo30';testgammar30'];
testerror30=[testden30';testporo30'];
fprintf(fid,'%6.4f %6.4f\n',testerror30);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror30);
fclose(fid);
 
fid = fopen('tinerror40.dat','w');
%testerror40=[testden40';testporo40';testgammar40'];
testerror40=[testden40';testporo40'];
fprintf(fid,'%6.4f %6.4f\n',testerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror40);
fclose(fid);
 
fid = fopen('tinerror50.dat','w');
%testerror50=[testden50';testporo50';testgammar50'];
testerror50=[testden50';testporo50'];
fprintf(fid,'%6.4f %6.4f\n',testerror50);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror50);
fclose(fid);
 
load tinerror10.dat;
load tinerror20.dat;
load tinerror30.dat;
load tinerror40.dat;
load tinerror50.dat
 
 
[e10t]= tramnmx(tinerror10',minp,maxp);
[p2e10t] = trapca(e10t,transMat);
a10t = gpfwd(net,p2e10t');
[a10t] = postmnmx(a10t',mint,maxt);
 
 
tat=tat';
tdev10=tat'-a10t;
tdev10=tdev10';
tdp10=tdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
tstmse10=mse(tdp10);
 
figure(45);
 
subplot(1,1,1),plot(tdp10,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
[e20t]= tramnmx(tinerror20',minp,maxp);
[p2e20t] = trapca(e20t,transMat);
a20t = gpfwd(net,p2e20t');
[a20t] = postmnmx(a20t',mint,maxt);
 
 
 
tdev20=tat'-a20t;
tdev20=tdev20';
tdp20=tdev20(:,1);
tstmse20=mse(tdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(46);
 
subplot(1,1,1),plot(tdp20,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
[e30t]= tramnmx(tinerror30',minp,maxp);
[p2e30t] = trapca(e30t,transMat);
a30t = gpfwd(net,p2e10t');
[a30t] = postmnmx(a30t',mint,maxt);
 
 
 
tdev30=tat'-a30t;
tdev30=tdev30';
tdp30=tdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
tstmse30=mse(tdp30);
 
figure(47);
 
subplot(1,1,1),plot(tdp30,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
[e40t]= tramnmx(tinerror40',minp,maxp);
[p2e40t] = trapca(e40t,transMat);
a40t = gpfwd(net,p2e40t');
[a40t] = postmnmx(a40t',mint,maxt);
 
tdev40=tat'-a40t;
tdev40=tdev40';
tdp40=tdev40(:,1);
tstmse40=mse(tdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(48);
 
subplot(1,1,1),plot(tdp40,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
[e50t]= tramnmx(tinerror50',minp,maxp);
[p2e50t] = trapca(e50t,transMat);
a50t = gpfwd(net,p2e50t');
[a50t] = postmnmx(a50t',mint,maxt);
 
 
 
tdev50=tat'-a50t;
tdev50=tdev50';
tdp50=tdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
tstmse50=mse(tdp50);
 
figure(50);
 
subplot(1,1,1),plot(tdp50,'ro-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',20);
ylabel('Error Dev.','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
 
figure(51)
nx=[10 20 30 40 50];
msetst=[tstmse10 tstmse20 tstmse30 tstmse40 tstmse50];
bar(nx,msetst);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Test Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Obs.','GPreg', 'GPreg+STD','GPreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
figure(52)
count=[msetr' mseval' msetst']
ybar=count(1:5,:);
bar(ybar,'LineWidth',3,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('MSE Analysis in Noisy Input Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
meancount=mean(count')
figure(53)
count=[mseval' msetst'];
ybar=count(1:5,:);
bar(ybar,'LineWidth',3,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('MSE Analysis in Noisy Input Data');
legend('Validation','Test','Location','northeast');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
plotregression(tr,a10tr')
xlabel('Target','FontName','Tahoma','Fontsize',20);ylabel('GPreg','FontName','Tahoma','Fontsize',20);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Training Data+10% Noise','Location','northeast');
 
 
plotregression(tav,a10v')
xlabel('Target','FontName','Tahoma','Fontsize',20);ylabel('GPreg','FontName','Tahoma','Fontsize',20);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
 
 
plotregression(tat,a10t')
xlabel('Target','FontName','Tahoma','Fontsize',20);ylabel('GPreg','FontName','Tahoma','Fontsize',20);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Test Data+10% Noise','Location','northeast');
%correlation analysis at different target noise
vinerror10=[ff,gg];
[e10]= tramnmx(vinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
a10 = gpfwd(net,p2e10');
[a10] = postmnmx(a10',mint,maxt);
 
noisetarget10tr=target10(:,iitr);
noisetarget20tr=target20(:,iitr);
noisetarget30tr=target30(:,iitr);
noisetarget40tr=target40(:,iitr);
noisetarget50tr=target50(:,iitr);
 
noisetarget10tv=target10(:,iival);
noisetarget20tv=target20(:,iival);
noisetarget30tv=target30(:,iival);
noisetarget40tv=target40(:,iival);
noisetarget50tv=target50(:,iival);
 
noisetarget10tst=target10(:,iitst);
noisetarget20tst=target20(:,iitst);
noisetarget30tst=target30(:,iitst);
noisetarget40tst=target40(:,iitst);
noisetarget50tst=target50(:,iitst);
 
noisea10tr=a10(:,iitr);
noisea10tv=a10(:,iival);
noisea10tst=a10(:,iitst);
 
%%%%%%%%%% R analysis
trcorrcoef10=corrcoef(noisetarget10tr,noisea10tr);
valcorrcoef10=corrcoef(noisetarget10tv,noisea10tv);
testcorrcoef10=corrcoef(noisetarget10tst,noisea10tst);
 
noisetrcod10=diag(trcorrcoef10,1).*diag(trcorrcoef10,1);
noisevalcod10=diag(valcorrcoef10,1).*diag(valcorrcoef10,1);
noisetestcod10=diag(testcorrcoef10,1).*diag(testcorrcoef10,1);
%%%%
trcorrcoef20=corrcoef(noisetarget20tr,noisea10tr);
valcorrcoef20=corrcoef(noisetarget20tv,noisea10tv);
testcorrcoef20=corrcoef(noisetarget20tst,noisea10tst);
 
noisetrcod20=diag(trcorrcoef20,1).*diag(trcorrcoef20,1);
noisevalcod20=diag(valcorrcoef20,1).*diag(valcorrcoef20,1);
noisetestcod20=diag(testcorrcoef20,1).*diag(testcorrcoef20,1);
%%%%%%%%%%
trcorrcoef30=corrcoef(noisetarget30tr,noisea10tr);
valcorrcoef30=corrcoef(noisetarget30tv,noisea10tv);
testcorrcoef30=corrcoef(noisetarget30tst,noisea10tst);
 
noisetrcod30=diag(trcorrcoef30,1).*diag(trcorrcoef30,1);
noisevalcod30=diag(valcorrcoef30,1).*diag(valcorrcoef30,1);
noisetestcod30=diag(testcorrcoef30,1).*diag(testcorrcoef30,1);
%%%%%%%%%%%%%%%%%%%5
trcorrcoef40=corrcoef(noisetarget40tr,noisea10tr);
valcorrcoef40=corrcoef(noisetarget40tv,noisea10tv);
testcorrcoef40=corrcoef(noisetarget40tst,noisea10tst);
 
noisetrcod40=diag(trcorrcoef40,1).*diag(trcorrcoef40,1);
noisevalcod40=diag(valcorrcoef40,1).*diag(valcorrcoef40,1);
noisetestcod40=diag(testcorrcoef40,1).*diag(testcorrcoef40,1);
%%%%%%
trcorrcoef50=corrcoef(noisetarget50tr,noisea10tr);
valcorrcoef50=corrcoef(noisetarget50tv,noisea10tv);
testcorrcoef50=corrcoef(noisetarget50tst,noisea10tst);
 
noisetrcod50=diag(trcorrcoef50,1).*diag(trcorrcoef50,1);
noisevalcod50=diag(valcorrcoef50,1).*diag(valcorrcoef50,1);
noisetestcod50=diag(testcorrcoef50,1).*diag(testcorrcoef50,1);
 
r2tr=[noisetrcod10 noisetrcod20 noisetrcod30 noisetrcod40 noisetrcod50];
r2tv=[noisevalcod10 noisevalcod20 noisevalcod30 noisevalcod40 noisevalcod50];
r2tst=[noisetestcod10 noisetestcod20 noisetestcod30 noisetestcod40 noisetestcod50];
 
r2tr=sqrt(r2tr);
r2tv=sqrt(r2tv);
r2tst=sqrt(r2tst);
 
figure(70)
countr2=[r2tr' r2tv' r2tst']
ybar2=countr2(1:5,:);
bar(ybar2,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('R ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('R Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2=mean(countr2')
%..................MSE analysis
trmsecorrcoef10=mse(noisetarget10tr-noisea10tr);
valmsecorrcoef10=mse(noisetarget10tv-noisea10tv);
testmsecorrcoef10=mse(noisetarget10tst-noisea10tst);
 
trmsecorrcoef20=mse(noisetarget20tr-noisea10tr);
valmsecorrcoef20=mse(noisetarget20tv-noisea10tv);
testmsecorrcoef20=mse(noisetarget20tst-noisea10tst);
 
trmsecorrcoef30=mse(noisetarget30tr-noisea10tr);
valmsecorrcoef30=mse(noisetarget30tv-noisea10tv);
testmsecorrcoef30=mse(noisetarget30tst-noisea10tst);
 
trmsecorrcoef40=mse(noisetarget40tr-noisea10tr);
valmsecorrcoef40=mse(noisetarget40tv-noisea10tv);
testmsecorrcoef40=mse(noisetarget40tst-noisea10tst);
 
trmsecorrcoef50=mse(noisetarget50tr-noisea10tr);
valmsecorrcoef50=mse(noisetarget50tv-noisea10tv);
testmsecorrcoef50=mse(noisetarget50tst-noisea10tst);
 
r2msetr=[trmsecorrcoef10 trmsecorrcoef20 trmsecorrcoef30 trmsecorrcoef40 trmsecorrcoef50];
r2msetv=[valmsecorrcoef10 valmsecorrcoef20 valmsecorrcoef30 valmsecorrcoef40 valmsecorrcoef50];
r2msetst=[testmsecorrcoef10 testmsecorrcoef20 testmsecorrcoef30 testmsecorrcoef40 testmsecorrcoef50];
 
figure(71)
countr2mse=[r2msetr' r2msetv' r2msetst']
ybar2mse=countr2mse(1:5,:);
bar(ybar2mse,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('MSE Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2mse=mean(countr2mse')
%Input noise analysis based on R analysis
%%%%%%%%%% R analysis
in10tr=a10tr;
in20tr=a20tr;
in30tr=a30tr;
in40tr=a40tr;
in50tr=a50tr;
 
in10tv=a10v;
in20tv=a20v;
in30tv=a30v;
in40tv=a40v;
in50tv=a50v;
 
in10tst=a10t;
in20tst=a20t;
in30tst=a30t;
in40tst=a40t;
in50tst=a50t;
 
out10tr=tr;
out20tr=tr;
out30tr=tr;
out40tr=tr;
out50tr=tr;
 
out10tv=tav;
out20tv=tav;
out30tv=tav;
out40tv=tav;
out50tv=tav;
 
out10tst=tat;
out20tst=tat;
out30tst=tat;
out40tst=tat;
out50tst=tat;
 
trcorrcoef10=corrcoef(in10tr,out10tr);
valcorrcoef10=corrcoef(in10tv,out10tv);
testcorrcoef10=corrcoef(in10tst,out10tst);
 
noisetrcod10=diag(trcorrcoef10,1).*diag(trcorrcoef10,1);
noisevalcod10=diag(valcorrcoef10,1).*diag(valcorrcoef10,1);
noisetestcod10=diag(testcorrcoef10,1).*diag(testcorrcoef10,1);
%%%%
trcorrcoef20=corrcoef(in20tr,out20tr);
valcorrcoef20=corrcoef(in20tv,out20tv);
testcorrcoef20=corrcoef(in20tst,out20tst);
 
noisetrcod20=diag(trcorrcoef20,1).*diag(trcorrcoef20,1);
noisevalcod20=diag(valcorrcoef20,1).*diag(valcorrcoef20,1);
noisetestcod20=diag(testcorrcoef20,1).*diag(testcorrcoef20,1);
%%%%%%%%%%
trcorrcoef30=corrcoef(in30tr,out30tr);
valcorrcoef30=corrcoef(in30tv,out30tv);
testcorrcoef30=corrcoef(in30tst,out30tst);
 
noisetrcod30=diag(trcorrcoef30,1).*diag(trcorrcoef30,1);
noisevalcod30=diag(valcorrcoef30,1).*diag(valcorrcoef30,1);
noisetestcod30=diag(testcorrcoef30,1).*diag(testcorrcoef30,1);
%%%%%%%%%%%%%%%%%%%5
trcorrcoef40=corrcoef(in40tr,out40tr);
valcorrcoef40=corrcoef(in40tv,out40tv);
testcorrcoef40=corrcoef(in40tst,out40tst);
 
noisetrcod40=diag(trcorrcoef40,1).*diag(trcorrcoef40,1);
noisevalcod40=diag(valcorrcoef40,1).*diag(valcorrcoef40,1);
noisetestcod40=diag(testcorrcoef40,1).*diag(testcorrcoef40,1);
%%%%%%
trcorrcoef50=corrcoef(in50tr,out50tr);
valcorrcoef50=corrcoef(in50tv,out50tv);
testcorrcoef50=corrcoef(in50tst,out50tst);
 
noisetrcod50=diag(trcorrcoef50,1).*diag(trcorrcoef50,1);
noisevalcod50=diag(valcorrcoef50,1).*diag(valcorrcoef50,1);
noisetestcod50=diag(testcorrcoef50,1).*diag(testcorrcoef50,1);
 
r2tr=[noisetrcod10 noisetrcod20 noisetrcod30 noisetrcod40 noisetrcod50];
r2tv=[noisevalcod10 noisevalcod20 noisevalcod30 noisevalcod40 noisevalcod50];
r2tst=[noisetestcod10 noisetestcod20 noisetestcod30 noisetestcod40 noisetestcod50];
 
r2tr=sqrt(r2tr)
r2tv=sqrt(r2tv)
r2tst=sqrt(r2tst)
 
figure(72)
countr2=[r2tr' r2tv' r2tst']
ybar2=countr2(1:5,:);
bar(ybar2,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('R ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('R Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2=mean(countr2')


% % ARD-BNN  for the MLP...............................................
% % pause;
% 
clc;
clear all;
close all;
 
%Training data;
load gravdeep.dat;
aa=gravdeep(:,1);%longitude
bb=gravdeep(:,2);%latitude
cc=gravdeep(:,3);%BG value
 
dd=gravdeep(:,4);%Altitude
ee=gravdeep(:,5);%deep basement
% %----------------------------------------------------
%Training data;
% load gravshalow.dat;%205 no.
% a=gravshalow(:,1);%longitude
% b=gravshalow(:,2);%latitude
% c=gravshalow(:,3);%BG value
% d=gravshalow(:,4);%altitude
% e=gravshalow(:,5);%shallow basement
%%*************************TRAINING NOISY DATA****************************************
% aa=aa';
% bb=bb';
% cc=cc';
 
totalden = aa;
totalporo = bb;
totalgammar = cc;
%........................................
     
%.........................................
%White Noise stability test
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
Aa1=0.99;%0.99
dsm(1)=a(1);
for i=2:nmax,
dsm(i)=Aa1.*dsm(i-1)+r(i);
end
dsm=dsm/max(abs(dsm));
dsm=dsm';
%r=r';
%den0=den.*r;
%den=a;
%totalden0=dsm;
totalden0=totalden.*dsm;
totalden10=totalden+0.001.*totalden0;
totalden20=totalden+0.002.*totalden0;
totalden30=totalden+0.003.*totalden0;
totalden40=totalden+0.004.*totalden0;
totalden50=totalden+0.005.*totalden0;
 
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
pr(1)=r(1);
Ba2=-0.0957;%-0.0957
for i=2:nmax,
pr(i)=Ba2.*pr(i-1)+r(i);
end
 
pr=pr/max(abs(pr));
pr=pr';
%poro0=poro.*r;
%poro=b;
%totalporo0=pr;
totalporo0=totalporo.*pr;
totalporo10=totalporo+0.001.*totalporo0;
totalporo20=totalporo+0.002.*totalporo0;
totalporo30=totalporo+0.003.*totalporo0;
totalporo40=totalporo+0.004.*totalporo0;
totalporo50=totalporo+0.005.*totalporo0;
 
randn('state',0);
nmax=211;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
gg(1)=r(1);
Ca3=0.0142;%0.0142
for i=2:nmax,
gg(i)=Ca3.*gg(i-1)+r(i);
end
gg=gg/max(abs(gg));
gg=gg';
%gammar0=gammar.*r;
%gammar=gammar';
%totalgammar0=gg;
totalgammar0=totalgammar.*gg;
totalgammar10=totalgammar+0.001.*totalgammar0;
totalgammar20=totalgammar+0.002.*totalgammar0;
totalgammar30=totalgammar+0.003.*totalgammar0;
totalgammar40=totalgammar+0.004.*totalgammar0;
totalgammar50=totalgammar+0.005.*totalgammar0;
 
fid = fopen('totalinerror10.dat','w');
%totalerror10=[totalden10';totalporo10';totalgammar10'];
totalerror10=[totalden10';totalporo10'];
fprintf(fid,'%6.4f %6.4f\n',totalerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror10);
fclose(fid);
 
fid = fopen('totalinerror20.dat','w');
%totalerror20=[totalden20';totalporo20';totalgammar20'];
totalerror20=[totalden20';totalporo20'];
fprintf(fid,'%6.4f %6.4f\n',totalerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror20);
fclose(fid);
 
fid = fopen('totalinerror30.dat','w');
%totalerror30=[totalden30';totalporo30';totalgammar30'];
totalerror30=[totalden30';totalporo30';];
fprintf(fid,'%6.4f %6.4f\n',totalerror30);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror30);
fclose(fid);
 
fid = fopen('totalinerror40.dat','w');
%totalerror40=[totalden40';totalporo40';totalgammar40'];
totalerror40=[totalden40';totalporo40'];
fprintf(fid,'%6.4f %6.4f\n',totalerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror40);
fclose(fid);
 
fid = fopen('totalinerror50.dat','w');
%totalerror50=[totalden50';totalporo50';totalgammar50'];
totalerror50=[totalden50';totalporo50'];
fprintf(fid,'%6.4f %6.4f\n',totalerror50);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror50);
fclose(fid);
 
load totalinerror10.dat
load totalinerror20.dat
load totalinerror30.dat
load totalinerror40.dat
load totalinerror50.dat
 
%Training data;
totalinerror=totalinerror30;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
%ccn=totalinerror(:,3);%BG value
 
%........................................................................
% Real data
load gravwhole.dat;
ff=gravwhole(:,1);%longitude
gg=gravwhole(:,2);%latitude
hh=gravwhole(:,3);%BG value
ii=gravwhole(:,4);%altitude
 
x=[aan,bbn];
t=ee;
 
x=x';
t=t';
 
%.....................................................noise adding......
 %ndata=211;
 noise = 0.15;          % Standard deviation of noise distribution.
 
 %*********************************************************
target=t;
rand('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
Ta1=0.99;%0.99
dm(1)=a(1);
for i=2:nmax,
dm(i)=Ta1.*dm(i-1)+r(i);
end
 
dm=dm./max(abs(dm));
%dm=dm';
%r=r';
%den0=den.*r;
%den=a;
%target0=dm';
target0=target.*dm;
target10=target+0.1.*target0;
target20=target+0.2.*target0;
target30=target+0.3.*target0;
target40=target+0.4.*target0;
target50=target+0.5.*target0;
t=target30;
 
 
%**********************************************************************
 
%***********************************************************************
[xn,minp,maxp,tn,mint,maxt] = premnmx(x,t);
[xtrans,transMat] = prepca(xn,0.0000000000000002);
[R,Q] = size(xtrans);
iitst = 2:4:Q;
iival = 4:4:Q;
iitr = [1:4:Q 3:4:Q];
valP = xtrans(:,iival); val.T = tn(:,iival);
testP = xtrans(:,iitst); test.T = tn(:,iitst);
xtr = xtrans(:,iitr); ttr = tn(:,iitr);
xr=xtr';
zr=ttr';
 
 
 
% % pause;
 
% Set up network parameters.
nin = 2;        % Number of inputs.
nhidden = 15;       % Number of hidden units.
nout = 1;       % Number of outputs.
alpha = 0.01;       % Initial prior hyperparameter. 
beta_init = 50; % Initial noise hyperparameter.
 
% Create and initialize network weight vector.
net = mlp(nin, nhidden, nout, 'linear', alpha, beta_init);
 
% Set up vector of options for the optimiser.
nouter = 5;         % Number of outer loops.
ninner = 5;         % Number of innter loops.
options = zeros(1,20);      % Default options vector.
options(1) = 1;         % This provides display of error values.
options(2) = 1.0e-7;        % Absolute precision for weights.
options(3) = 1.0e-7;        % Precision for objective function.
options(14) = 500;      % Number of training cycles in inner loop. 
 
% Train using scaled conjugate gradients, re-estimating alpha and beta.
for k = 1:nouter
  net = netopt(net, options, xr, zr, 'scg');
  [net, gamma] = evidence(net, xr, zr, ninner);
  fprintf(1, '\nRe-estimation cycle %d:\n', k);
  fprintf(1, '  alpha =  %8.5f\n', net.alpha);
  fprintf(1, '  beta  =  %8.5f\n', net.beta);
  fprintf(1, '  gamma =  %8.5f\n\n', gamma);
  disp(' ')
  disp('Press any key to continue.')
%   pause;
end
 
 fprintf(1, 'true beta: %f\n', 1/(noise*noise));
 
%...............................
[ymaint, mainsigt] = netevfwd(mlppak(net), net, x', t', x');
sigmaint = sqrt(mainsigt);
sigmaint=abs(sigmaint);
sigtr=sigmaint([1:4:Q 3:4:Q]);
msigtr=mean(sigtr);
sigval=sigmaint(2:4:Q);
msigval=mean(sigval);
sigtst=sigmaint(4:4:Q);
msigtst=mean(sigtst);
 
 
sqxr=sigtr;
[an] = mlpfwd(net,xr);
[y] = postmnmx(an',mint,maxt);
sqtr=2.00*(sqxr');
tr=t(:,iitr);
 
figure(26);
 
xis=1:106;
E1=sqtr(:,1).*ones(size(y));
%subplot(3,1,1),
errorbar(xis,y,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tr,'ro-',xis, y,'o-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -6 6]);
legend('Error-bar','Target','ARD-BNNreg.');
xlabel('No. of Data','Fontsize',24);
ylabel('Sediment Depth(Km)','Fontsize',24);
title('Training Interval')
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
%........................................
figure;
an = mlpfwd(net,xr);
[y] = postmnmx(an',mint,maxt);
ee=ee';
tr=ee(:,iitr);
plot(tr,'-r','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
plot(y,'--b>','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
 plot(y + msigtr, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
 plot(y - msigtr, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
xlabel('Number of Training Samples','FontName','Tahoma','Fontsize',24);ylabel('Shallow Interface Depth(Km)','FontName','Tahoma','Fontsize',24);
title('Shallow Interface Depth Estimation:Training Interval','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
legend('Obs.','ARD-BNNreg', 'ARD-BNNreg+STD','ARD-BNNreg-STD');
 
plotregression(tr,y)
xlabel('Obs.','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Basement Depth Prediction:Training Interval','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
 

%............................................
sqvalx=sigval;
[av] = mlpfwd(net,valP');
[avy] = postmnmx(av',mint,maxt);
sqval=2.0*(sqvalx');
[av] = postmnmx(av',mint,maxt);
%ee=ee';
tav=ee(:,iival);
 
figure(27);
xis=1:52;
E2=sqtr(:,1).*ones(size(av));
 
%subplot(3,1,1),
errorbar(xis,avy,E2,'c');grid on;axis([0 52 0 6]);hold on
%title('Errobar plot');hold on
plot(xis,tav,'ro-',xis, avy,'o-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 52 0 6]);
legend('Error-bar','Target','ARD-BNNreg.');
%plot(tat);hold on;plot(ate);
xlabel('No. of Data','Fontsize',24);
ylabel('Sediment Depth(Km)','Fontsize',24);
title('Validation Interval')
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
  
figure;
av = mlpfwd(net,valP');
[av] = postmnmx(av',mint,maxt);
%tav=ee(:,iival);
plot(tav,'-r','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
plot(av,'--b>','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
 plot(av + msigval, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
 plot(av - msigval, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
xlabel('Number of Validation Samples','FontName','Tahoma','Fontsize',24);ylabel('Shallow Interface Depth(Km) ','FontName','Tahoma','Fontsize',24);
title('Shallow Interface Depth Estimation:Validation Interval','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
legend('Obs.','ARD-BNNreg', 'ARD-BNNreg+STD','ARD-BNNreg-STD');
 
 
plotregression(tav,av)
xlabel('Obs.','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Basement Depth Prediction:Training Interval','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
 
sqtestx=sigtst;
[ate] = mlpfwd(net,testP');
[atey] = postmnmx(ate',mint,maxt);
sqtest=2.0*(sqtestx);
%tat=t(:,iitst);
tat=ee(:,iitst);
 
figure(28);
xis=1:53;
E3=sqtr(:,1).*ones(size(atey));
 
%subplot(3,1,1),
errorbar(xis,atey,E3,'c');grid on;axis([0 53 0 6]);hold on
%title('Errobar plot');hold on
plot(xis,tat,'ro-',xis, atey,'o-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 53 0 6]);
legend('Error-bar','Target','ARD-BNNreg.');
%plot(tat);hold on;plot(ate);
xlabel('No. of Data','Fontsize',24);
ylabel('Sediment Depth(Km)','Fontsize',24);
title('Test Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
figure
ate = mlpfwd(net,testP');
[ate] = postmnmx(ate',mint,maxt);
%tat=ee(:,iitst);
plot(tat,'-r','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on
plot(ate,'--b>','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6); hold on;
 plot(ate + msigtst, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);hold on;
 plot(ate - msigtst, '-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
xlabel('Number of test Samples','FontName','Tahoma','Fontsize',24);ylabel('Shallow Interface Depth(Km) ','FontName','Tahoma','Fontsize',24);
title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
legend('Obs.','ARD-BNNreg', 'ARD-BNNreg+STD','ARD-BNNreg-STD');
 
 
plotregression(tat,ate)
xlabel('Obs.','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Basement Depth Prediction:Training Interval','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
 
figure
% Real data analysis;
vinerror10=[ff,gg];
[e10]= tramnmx(vinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
a10 = mlpfwd(net,p2e10');
[a10] = postmnmx(a10',mint,maxt);
%[a10]=-a10;
 
fdata=488;
gdata=488;
fmax=max(ff);
gmax=max(gg);
fmin=min(ff);
gmin=min(gg);
fv = linspace(fmin, fmax,fdata);
gv = linspace(gmin, gmax,gdata);
[fi,gi] = meshgrid(fv,gv);
c1i = griddata(ff,gg,a10,fi,gi,'cubic');
 
meshz(fi,gi,c1i);hold on;
colorbar;
plot3(ff,gg,a10,'v','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
%view(10,90);
% view(0,90);
%shading interp;
grid on;
box on;
%axis tight;
%axis auto;
%axis square
colorbar;
set(gca,'LineWidth',2,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('Longitude(in degree)','Fontsize',24);
ylabel('Latitude(in degree)','Fontsize',24);
zlabel('Deeper Interface Depth(Km)','Fontsize',24);set(gca,'Zdir','reverse');
title('Deeper Interface Depth Estimation','Fontsize',24);
%.........................................................................
figure;
[ymain, mainsig] = netevfwd(mlppak(net), net, x', t', vinerror10);
sigmain = sqrt(mainsig);
sigmain=abs(sigmain);
load datav2.dat;
sigmain1=datav2(:,1);
 
c2i = griddata(ff,gg,sigmain1,fi,gi,'cubic');
 
meshz(fi,gi,c2i);hold on;%meshc(fi,gi,c1i);hold on;
%surf(fi,gi,c1i);hold on;
%surf(ai,bi,c1i);hold on;
%[cout,H,cf]=contourf(fi,gi,c1i);
%f(cout,H);hold on
plot3(ff,gg,sigmain1,'v','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
%view(10,90);
% view(0,90);
%shading interp;
grid on;
box on;
%axis tight;
%axis auto;
%axis square
colorbar;
set(gca,'LineWidth',2,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('Longitude(in degree)','Fontsize',24);
ylabel('Latitude(in degree)','Fontsize',24);
zlabel('STD (Km)','Fontsize',24);%set(gca,'Zdir','reverse');
title('Uncertainty Estimation in Deeper Interface Depth Prediction ','Fontsize',24);
%vt = t(:,iival);
 
disp('Press any key to end.')

%Statistical %Training....................................................
%.........................................................................
ntotal=211;
ntrain=106;
nval=52;
ntest=53;
 
x1train=tr(1:ntrain);
x1val=tav(1:nval);
x1test=tat(1:ntest);
 
x2train=y(1:ntrain);
x2val=av(1:nval);
x2test=ate(1:ntest);
 
disp('Statistics of x1/obs......................................');
x1avgtrain=mean(x1train)
x1avgnval=mean(x1val)
x1avgtest=mean(x1test)
 
x1stdtrain=std(x1train)
x1stdval=std(x1val)
x1stdtest=std(x1test)
 
x1trmedian=median(x1train)
x1valmedian=median(x1val)
x1testmedian=median(x1test)
 
x1trmode=mode(x1train)
x1valmode=mode(x1val)
x1testmode=mode(x1test)
 
x1trskewness=skewness(x1train)
x1valskewness=skewness(x1val)
x1testskewness=skewness(x1test)
 
x1trkurtosis=kurtosis(x1train)
x1valkurtosis=kurtosis(x1val)
x1testkurtosis=kurtosis(x1test)
 
disp('Statistics of x2/ARD-BNNreg......................................')
x2avgtrain=mean(x2train)
x2avgnval=mean(x2val)
x2avgtest=mean(x2test)
 
x2stdtrain=std(x2train)
x2stdval=std(x2val)
x2stdtest=std(x2test)
 
x2trmedian=median(x2train);
x2valmedian=median(x2val)
x2testmedian=median(x2test)
 
x2trmode=mode(x2train)
x2valmode=mode(x2val)
x2testmode=mode(x2test)
 
x2trskewness=skewness(x2train)
x2valskewness=skewness(x2val)
x2testskewness=skewness(x2test)
 
x2trkurtosis=kurtosis(x2train)
x2valkurtosis=kurtosis(x2val)
x2testkurtosis=kurtosis(x2test)
 
x2errtrain=mse(x2train-x1train)
x2errval=mse(x2val-x1val)
x2errtest=mse(x2test-x1test)
 
x2retrain=1.0-sum((x1train-x2train).^2)./sum(x1train.^2)
x2reval=1.0-sum((x1val-x2val).^2)./sum(x1val.^2)
x2retest=1.0-sum((x1test-x2test).^2)./sum(x1test.^2)
 
x2dtrain=1.0-sum((x1train-x2train).^2)./sum(abs(x2train-mean(x1train))+abs(x1train-mean(x1train)).^2)
x2dval=1.0-sum((x1val-x2val).^2)./sum(abs(x2val-mean(x1val))+abs(x1val-mean(x1val)).^2)
x2dtest=1.0-sum((x1test-x2test).^2)./sum(abs(x2test-mean(x1test))+abs(x1test-mean(x1test)).^2)
 
trcorrcoef=corrcoef(x1train,x2train)
valcorrcoef=corrcoef(x1val,x2val)
testcorrcoef=corrcoef(x1test,x2test)
 
trcod=diag(trcorrcoef,1).*diag(trcorrcoef,1)
valcod=diag(valcorrcoef,1).*diag(valcorrcoef,1)
testcod=diag(testcorrcoef,1).*diag(testcorrcoef,1)
 
%***************Training Error Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error analysis**********************************************************************
aa=aa';
bb=bb';
cc=cc';
 
inputden = aa(:,iitr);
inputporo = bb(:,iitr);
inputgammar = cc(:,iitr);
 
inputden=inputden';
inputporo=inputporo';
inputgammar=inputgammar';
 
rand('state',0);
nmax=106;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
a1=0.94;%0.94
dtr(1)=a(1);
for i=2:nmax,
dtr(i)=a1.*dtr(i-1)+r(i);
end
 
dtr=dtr/max(abs(dtr));
dtr=dtr';
%r=r';
%den0=den.*r;
%den=a;
%inputden0=dtr;
inputden0=inputden.*dtr;
inputden10=inputden+0.001.*inputden0;
inputden20=inputden+0.002.*inputden0;
inputden30=inputden+0.003.*inputden0;
inputden40=inputden+0.004.*inputden0;
inputden50=inputden+0.005.*inputden0;
 
randn('state',0);
nmax=106;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
ptr(1)=r(1);
a2=-0.107;%-0.107
for i=2:nmax,
ptr(i)=a2.*ptr(i-1)+r(i);
end
ptr=ptr';
ptr=ptr/max(abs(ptr));
%poro0=poro.*r;
%poro=b;
%inputporo0=ptr;
inputporo0=inputporo.*ptr;
inputporo10=inputporo+0.001.*inputporo0;
inputporo20=inputporo+0.002.*inputporo0;
inputporo30=inputporo+0.003.*inputporo0;
inputporo40=inputporo+0.004.*inputporo0;
inputporo50=inputporo+0.005.*inputporo0;
 
randn('state',0);
nmax=106;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
ggtr(1)=r(1);
a3=0.122;%0.122
for i=2:nmax,
ggtr(i)=a3.*ggtr(i-1)+r(i);
end
 
ggtr=ggtr/max(abs(ggtr));
ggtr=ggtr';
%gammar0=gammar.*r;
%gammar=gammar';
%inputgammar0=ggtr;
inputgammar0=inputgammar.*ggtr;
inputgammar10=inputgammar+0.001.*inputgammar0;
inputgammar20=inputgammar+0.002.*inputgammar0;
inputgammar30=inputgammar+0.003.*inputgammar0;
inputgammar40=inputgammar+0.004.*inputgammar0;
inputgammar50=inputgammar+0.005.*inputgammar0;
 
fid = fopen('sinerror10.dat','w');
%inputerror10=[inputden10';inputporo10';inputgammar10'];
inputerror10=[inputden10';inputporo10'];
fprintf(fid,'%6.4f %6.4f\n',inputerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror10);
fclose(fid);
 
fid = fopen('sinerror20.dat','w');
%inputerror20=[inputden20';inputporo20';inputgammar20'];
inputerror20=[inputden20';inputporo20'];
fprintf(fid,'%6.4f %6.4f\n',inputerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror20);
fclose(fid);
 
fid = fopen('sinerror30.dat','w');
%inputerror30=[inputden30';inputporo30';inputgammar30'];
inputerror30=[inputden30';inputporo30'];
fprintf(fid,'%6.4f %6.4f\n',inputerror30);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror30);
fclose(fid);
 
fid = fopen('sinerror40.dat','w');
%inputerror40=[inputden40';inputporo40';inputgammar40'];
inputerror40=[inputden40';inputporo40'];
fprintf(fid,'%6.4f %6.4f\n',inputerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror40);
fclose(fid);
 
fid = fopen('sinerror50.dat','w');
%inputerror50=[inputden50';inputporo50';inputgammar50'];
inputerror50=[inputden50';inputporo50'];
fprintf(fid,'%6.4f %6.4f\n',inputerror50);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror50);
fclose(fid);
 
load sinerror10.dat;
load sinerror20.dat;
load sinerror30.dat;
load sinerror40.dat;
load sinerror50.dat
 
 
[e10]= tramnmx(sinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
a10 = mlpfwd(net,p2e10');
[a10tr] = postmnmx(a10',mint,maxt);
 
  
%t=t(:,iitr);
%tr=t(:,iitr);
trdev10=tr-a10tr;
trdev10=trdev10';
trdp10=trdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
trmse10=mse(trdp10);
 
figure(29);
 
subplot(1,1,1),plot(trdp10,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e20tr]= tramnmx(sinerror20',minp,maxp);
[p2e20tr] = trapca(e20tr,transMat);
a20tr = mlpfwd(net,p2e20tr');
[a20tr] = postmnmx(a20tr',mint,maxt);
 
 
 
trdev20=tr-a20tr;
trdev20=trdev20';
trdp20=trdev20(:,1);
trmse20=mse(trdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(30);
 
subplot(1,1,1),plot(trdp20,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e30tr]= tramnmx(sinerror30',minp,maxp);
[p2e30tr] = trapca(e30tr,transMat);
a30tr = mlpfwd(net,p2e30tr');
[a30tr] = postmnmx(a30tr',mint,maxt);
 
trdev30=tr-a30tr;
trdev30=trdev30';
trdp30=trdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
trmse30=mse(trdp30);
 
figure(31);
 
subplot(1,1,1),plot(trdp30,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e40tr]= tramnmx(sinerror40',minp,maxp);
[p2e40tr] = trapca(e40tr,transMat);
a40tr = mlpfwd(net,p2e40tr');
[a40tr] = postmnmx(a40tr',mint,maxt);
 
trdev40=tr-a40tr;
trdev40=trdev40';
trdp40=trdev40(:,1);
trmse40=mse(trdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(32);
 
subplot(1,1,1),plot(trdp10,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e50tr]= tramnmx(sinerror50',minp,maxp);
[p2e50tr] = trapca(e50tr,transMat);
a50tr = mlpfwd(net,p2e50tr');
[a50tr] = postmnmx(a50tr',mint,maxt);
 

trdev50=tr-a50tr;
trdev50=trdev50';
trdp50=trdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
trmse50=mse(trdp50);
 
figure(33);
subplot(1,1,1),plot(trdp50,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
%error diagram
figure(34)
 
nx=[10 20 30 40 50];
msetr=[trmse10 trmse20 trmse30 trmse40 trmse50];
bar(nx,msetr);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('Training Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%******************Validation Error ANALYSIS****************************
 
 
vdenv = aa(:,iival);
vporov = bb(:,iival);
vgammarv = cc(:,iival);
 
valden=vdenv';
valporo=vporov';
valgammar=vgammarv';
 
randn('state',0);
nmax=52;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
dv(1)=r(1);
A1=0.99;%0.99
for i=2:nmax,
dv(i)=A1.*dv(i-1)+r(i);
end
 
dv=dv/max(abs(dv));
dv=dv';
%r=r';
%den0=den.*r;
%valden0=dv;
valden0=valden.*dv;
valden10=valden+0.001.*valden0;
valden20=valden+0.002.*valden0;
valden30=valden+0.003.*valden0;
valden40=valden+0.004.*valden0;
valden50=valden+0.005.*valden0;
 
randn('state',0);
nmax=52;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
pv(1)=r(1);
A2=-0.0495;%-0.0495
for i=2:nmax,
pv(i)=A2.*pv(i-1)+r(i);
end
pv=pv';
pv=pv/max(abs(pv));
%poro0=poro.*r;
%poro=b;
%valporo0=pv;
valporo0=valporo.*pv;
valporo10=valporo+0.001.*valporo0;
valporo20=valporo+0.002.*valporo0;
valporo30=valporo+0.003.*valporo0;
valporo40=valporo+0.004.*valporo0;
valporo50=valporo+0.005.*valporo0;
 
randn('state',0);
nmax=52;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
gva(1)=r(1);
A3=0.1408;%0.1408
for i=2:nmax,
gva(i)=A3.*gva(i-1)+r(i);
end
 
gva=gva/max(abs(gva));
gva=gva';
valgammar0=valgammar.*gva;
valgammar10=valgammar+0.001.*valgammar0;
valgammar20=valgammar+0.002.*valgammar0;
valgammar30=valgammar+0.003.*valgammar0;
valgammar40=valgammar+0.004.*valgammar0;
valgammar50=valgammar+0.005.*valgammar0;
 
fid = fopen('vinerror10.dat','w');
%valerror10=[valden10';valporo10';valgammar10'];
valerror10=[valden10';valporo10'];
fprintf(fid,'%6.4f %6.4f\n',valerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror10);
fclose(fid);
 
fid = fopen('vinerror20.dat','w');
%valerror20=[valden20';valporo20';valgammar20'];
valerror20=[valden20';valporo20'];
fprintf(fid,'%6.4f %6.4f\n',valerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror20);
fclose(fid);
 
fid = fopen('vinerror30.dat','w');
%valerror30=[valden30';valporo30';valgammar30'];
valerror30=[valden30';valporo30'];
fprintf(fid,'%6.4f %6.4f\n',valerror30);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror30);
fclose(fid);
 
fid = fopen('vinerror40.dat','w');
%valerror40=[valden40';valporo40';valgammar40'];
valerror40=[valden40';valporo40'];
fprintf(fid,'%6.4f %6.4f\n',valerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror40);
fclose(fid);
 
fid = fopen('vinerror50.dat','w');
%valerror50=[valden50';valporo50';valgammar50'];
valerror50=[valden50';valporo50'];
fprintf(fid,'%6.4f %6.4f\n',valerror50);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror50);
fclose(fid);
 
load vinerror10.dat;
load vinerror20.dat;
load vinerror30.dat;
load vinerror40.dat;
load vinerror50.dat
 
 
[e10v]= tramnmx(vinerror10',minp,maxp);
[p2e10v] = trapca(e10v,transMat);
a10v = mlpfwd(net,p2e10v');
[a10v] = postmnmx(a10v',mint,maxt);
 
tav=tav';
vdev10=tav'-a10v;
vdev10=vdev10';
vdp10=vdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
valmse10=mse(vdp10);
 
figure(35);
 
subplot(1,1,1),plot(vdp10,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e20v]= tramnmx(vinerror20',minp,maxp);
[p2e20v] = trapca(e20v,transMat);
a20v = mlpfwd(net,p2e20v');
[a20v] = postmnmx(a20v',mint,maxt);
 
vdev20=tav'-a20v;
vdev20=vdev20';
vdp20=vdev20(:,1);
valmse20=mse(vdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(36);
 
subplot(1,1,1),plot(vdp20,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e30v]= tramnmx(vinerror30',minp,maxp);
[p2e30v] = trapca(e30v,transMat);
a30v = mlpfwd(net,p2e30v');
[a30v] = postmnmx(a30v',mint,maxt);
 
 
 
vdev30=tav'-a30v;
vdev30=vdev30';
vdp30=vdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
valmse30=mse(vdp30);
 
figure(37);
 
subplot(1,1,1),plot(vdp30,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e40v]= tramnmx(vinerror40',minp,maxp);
[p2e40v] = trapca(e40v,transMat);
a40v = mlpfwd(net,p2e40v');
[a40v] = postmnmx(a40v',mint,maxt);
 
vdev40=tav'-a40v;
vdev40=vdev40';
vdp40=vdev40(:,1);
valmse40=mse(vdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(38);
 
subplot(1,1,1),plot(vdp40,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e50v]= tramnmx(vinerror50',minp,maxp);
[p2e50v] = trapca(e50v,transMat);
a50v = mlpfwd(net,p2e50v');
[a50v] = postmnmx(a50v',mint,maxt);
 
 
 
vdev50=tav'-a50v;
vdev50=vdev50';
vdp50=vdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
valmse50=mse(vdp50);
 
figure(40);
subplot(1,1,1),plot(vdp50,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
%error diagram
figure(41)
nx=[10 20 30 40 50];
mseval=[valmse10 valmse20 valmse30 valmse40 valmse50];
bar(nx,mseval);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);

title('Validation Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%**********************Test Error Analysis*****************************
%Error Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error analysis
dent = aa(:,iitst);
porot = bb(:,iitst);
gammart = cc(:,iitst);
 
testden=dent';
testporo=porot';
testgammar=gammart';
 
randn('state',0);
nmax=53;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
dtst(1)=r(1);
B1=0.99;%0.99
 
for i=2:nmax,
dtst(i)=B1.*dtst(i-1)+r(i);
end
 
dtst=dtst/max(abs(dtst));
dtst=dtst';
%r=r';
%den0=den.*r;
%den=a;
%testden0=dtst;
testden0=testden.*dtst;
testden10=testden+0.001.*testden0;
testden20=testden+0.002.*testden0;
testden30=testden+0.003.*testden0;
testden40=testden+0.004.*testden0;
testden50=testden+0.005.*testden0;
 
randn('state',0);
nmax=53;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
pt(1)=r(1);
B2=-0.2317;%-0.2317
for i=2:nmax,
pt(i)=B2.*pt(i-1)+r(i);
end
 
pt=pt/max(abs(pt));
pt=pt';
%poro0=poro.*r;
%poro=b;
%testporo0=pt;
testporo0=testporo.*pt;
testporo10=testporo+0.001.*testporo0;
testporo20=testporo+0.002.*testporo0;
testporo30=testporo+0.003.*testporo0;
testporo40=testporo+0.004.*testporo0;
testporo50=testporo+0.005.*testporo0;
 
randn('state',0);
nmax=53;
% White time series
r=randn(1,nmax);
r=r/max(abs(r));
gt(1)=r(1);
B3=-0.036317;%-0.036317
for i=2:nmax,
gt(i)=B3.*gt(i-1)+r(i);
end
 
gt=gt/max(abs(gt));
gt=gt';
%gammar0=gammar.*r;
%gammar=gammar';
%testgammar0=gt;
testgammar0=testgammar.*gt;
testgammar10=testgammar+0.001.*testgammar0;
testgammar20=testgammar+0.002.*testgammar0;
testgammar30=testgammar+0.003.*testgammar0;
testgammar40=testgammar+0.004.*testgammar0;
testgammar50=testgammar+0.005.*testgammar0;
 
fid = fopen('tinerror10.dat','w');
%testerror10=[testden10';testporo10';testgammar10'];
testerror10=[testden10';testporo10'];
fprintf(fid,'%6.4f %6.4f\n',testerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror10);
fclose(fid);
 
fid = fopen('tinerror20.dat','w');
%testerror20=[testden20';testporo20';testgammar20'];
testerror20=[testden20';testporo20'];
fprintf(fid,'%6.4f %6.4f\n',testerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror20);
fclose(fid);
 
fid = fopen('tinerror30.dat','w');
%testerror30=[testden30';testporo30';testgammar30'];
testerror30=[testden30';testporo30'];
fprintf(fid,'%6.4f %6.4f\n',testerror30);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror30);
fclose(fid);
 
fid = fopen('tinerror40.dat','w');
%testerror40=[testden40';testporo40';testgammar40'];
testerror40=[testden40';testporo40'];
fprintf(fid,'%6.4f %6.4f\n',testerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror40);
fclose(fid);
 
fid = fopen('tinerror50.dat','w');
%testerror50=[testden50';testporo50';testgammar50'];
testerror50=[testden50';testporo50'];
fprintf(fid,'%6.4f %6.4f\n',testerror50);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror50);
fclose(fid);
 
load tinerror10.dat;
load tinerror20.dat;
load tinerror30.dat;
load tinerror40.dat;
load tinerror50.dat
 
 
[e10t]= tramnmx(tinerror10',minp,maxp);
[p2e10t] = trapca(e10t,transMat);
a10t = mlpfwd(net,p2e10t');
[a10t] = postmnmx(a10t',mint,maxt);
 
 
tat=tat';
tdev10=tat'-a10t;
tdev10=tdev10';
tdp10=tdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
tstmse10=mse(tdp10);
 
figure(45);
 
subplot(1,1,1),plot(tdp10,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e20t]= tramnmx(tinerror20',minp,maxp);
[p2e20t] = trapca(e20t,transMat);
a20t = mlpfwd(net,p2e20t');
[a20t] = postmnmx(a20t',mint,maxt);
 
 
 
tdev20=tat'-a20t;
tdev20=tdev20';
tdp20=tdev20(:,1);
tstmse20=mse(tdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(46);
 
subplot(1,1,1),plot(tdp20,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
[e30t]= tramnmx(tinerror30',minp,maxp);
[p2e30t] = trapca(e30t,transMat);
a30t = mlpfwd(net,p2e10t');
[a30t] = postmnmx(a30t',mint,maxt);
 
 
 
tdev30=tat'-a30t;
tdev30=tdev30';
tdp30=tdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
tstmse30=mse(tdp30);
 
figure(47);
 
subplot(1,1,1),plot(tdp30,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
[e40t]= tramnmx(tinerror40',minp,maxp);
[p2e40t] = trapca(e40t,transMat);
a40t = mlpfwd(net,p2e40t');
[a40t] = postmnmx(a40t',mint,maxt);
 
tdev40=tat'-a40t;
tdev40=tdev40';
tdp40=tdev40(:,1);
tstmse40=mse(tdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(48);
 
subplot(1,1,1),plot(tdp40,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
[e50t]= tramnmx(tinerror50',minp,maxp);
[p2e50t] = trapca(e50t,transMat);
a50t = mlpfwd(net,p2e50t');
[a50t] = postmnmx(a50t',mint,maxt);
 
 
 
tdev50=tat'-a50t;
tdev50=tdev50';
tdp50=tdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
tstmse50=mse(tdp50);
 
figure(50);
 
subplot(1,1,1),plot(tdp50,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
 
figure(51)
nx=[10 20 30 40 50];
msetst=[tstmse10 tstmse20 tstmse30 tstmse40 tstmse50];
bar(nx,msetst);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('Test Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Obs.','ARD-BNNreg', 'ARD-BNNreg+STD','ARD-BNNreg-STD');
%set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
figure(52)
count=[msetr' mseval' msetst']
ybar=count(1:5,:);
bar(ybar,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('MSE Analysis in Noisy Input Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancount=mean(count')
 
figure(53)
count=[mseval' msetst'];
ybar=count(1:5,:);
bar(ybar,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('MSE Analysis in Noisy Input Data');
legend('Validation','Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
figure(54)
plotregression(tr,a10tr')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Training Data+10% Noise','Location','northeast');
 
figure(55)
plotregression(tav,a10v')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
figure(56)
plotregression(tav,a20v')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
figure(57)
plotregression(tav,a30v')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
 
figure(58)
plotregression(tav,a40v')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
 
figure(59)
plotregression(tav,a50v')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
 
figure(60)
plotregression(tat,a10t')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
figure(61)
plotregression(tat,a20t')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
figure(62)
plotregression(tat,a30t')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
figure(63)
plotregression(tat,a40t')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
figure(64)
plotregression(tat,a50t')
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ARD-BNNreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',20);axis tight
%legend('Test Data+10% Noise','Location','northeast');
 
%correlation analysis at different target noise
vinerror10=[ff,gg];
[e10]= tramnmx(vinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
a10 = mlpfwd(net,p2e10');
[a10] = postmnmx(a10',mint,maxt);
 
noisetarget10tr=target10(:,iitr);
noisetarget20tr=target20(:,iitr);
noisetarget30tr=target30(:,iitr);
noisetarget40tr=target40(:,iitr);
noisetarget50tr=target50(:,iitr);
 
noisetarget10tv=target10(:,iival);
noisetarget20tv=target20(:,iival);
noisetarget30tv=target30(:,iival);
noisetarget40tv=target40(:,iival);
noisetarget50tv=target50(:,iival);
 
noisetarget10tst=target10(:,iitst);
noisetarget20tst=target20(:,iitst);
noisetarget30tst=target30(:,iitst);
noisetarget40tst=target40(:,iitst);
noisetarget50tst=target50(:,iitst);
 
noisea10tr=a10(:,iitr);
noisea10tv=a10(:,iival);
noisea10tst=a10(:,iitst);
 
%%%%%%%%%% R analysis
trcorrcoef10=corrcoef(noisetarget10tr,noisea10tr)
valcorrcoef10=corrcoef(noisetarget10tv,noisea10tv)
testcorrcoef10=corrcoef(noisetarget10tst,noisea10tst)
 
noisetrcod10=diag(trcorrcoef10,1).*diag(trcorrcoef10,1)
noisevalcod10=diag(valcorrcoef10,1).*diag(valcorrcoef10,1)
noisetestcod10=diag(testcorrcoef10,1).*diag(testcorrcoef10,1)
%%%%
trcorrcoef20=corrcoef(noisetarget20tr,noisea10tr)
valcorrcoef20=corrcoef(noisetarget20tv,noisea10tv)
testcorrcoef20=corrcoef(noisetarget20tst,noisea10tst)
 
noisetrcod20=diag(trcorrcoef20,1).*diag(trcorrcoef20,1)
noisevalcod20=diag(valcorrcoef20,1).*diag(valcorrcoef20,1)
noisetestcod20=diag(testcorrcoef20,1).*diag(testcorrcoef20,1)
%%%%%%%%%%
trcorrcoef30=corrcoef(noisetarget30tr,noisea10tr)
valcorrcoef30=corrcoef(noisetarget30tv,noisea10tv)
testcorrcoef30=corrcoef(noisetarget30tst,noisea10tst)
 
noisetrcod30=diag(trcorrcoef30,1).*diag(trcorrcoef30,1)
noisevalcod30=diag(valcorrcoef30,1).*diag(valcorrcoef30,1)
noisetestcod30=diag(testcorrcoef30,1).*diag(testcorrcoef30,1)
%%%%%%%%%%%%%%%%%%%5
trcorrcoef40=corrcoef(noisetarget40tr,noisea10tr)
valcorrcoef40=corrcoef(noisetarget40tv,noisea10tv)
testcorrcoef40=corrcoef(noisetarget40tst,noisea10tst)
 
noisetrcod40=diag(trcorrcoef40,1).*diag(trcorrcoef40,1)
noisevalcod40=diag(valcorrcoef40,1).*diag(valcorrcoef40,1)
noisetestcod40=diag(testcorrcoef40,1).*diag(testcorrcoef40,1)
%%%%%%
trcorrcoef50=corrcoef(noisetarget50tr,noisea10tr)
valcorrcoef50=corrcoef(noisetarget50tv,noisea10tv)
testcorrcoef50=corrcoef(noisetarget50tst,noisea10tst)
 
noisetrcod50=diag(trcorrcoef50,1).*diag(trcorrcoef50,1)
noisevalcod50=diag(valcorrcoef50,1).*diag(valcorrcoef50,1)
noisetestcod50=diag(testcorrcoef50,1).*diag(testcorrcoef50,1)
 
r2tr=[noisetrcod10 noisetrcod20 noisetrcod30 noisetrcod40 noisetrcod50];
r2tv=[noisevalcod10 noisevalcod20 noisevalcod30 noisevalcod40 noisevalcod50];
r2tst=[noisetestcod10 noisetestcod20 noisetestcod30 noisetestcod40 noisetestcod50];
 
r2tr=sqrt(r2tr)
r2tv=sqrt(r2tv)
r2tst=sqrt(r2tst)
 
figure(70)
countr2=[r2tr' r2tv' r2tst']
ybar2=countr2(1:5,:);
bar(ybar2,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('R ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('R Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2=mean(countr2')
%..................MSE analysis
trmsecorrcoef10=mse(noisetarget10tr-noisea10tr)
valmsecorrcoef10=mse(noisetarget10tv-noisea10tv)
testmsecorrcoef10=mse(noisetarget10tst-noisea10tst)
 
trmsecorrcoef20=mse(noisetarget20tr-noisea10tr)
valmsecorrcoef20=mse(noisetarget20tv-noisea10tv)
testmsecorrcoef20=mse(noisetarget20tst-noisea10tst)
 
trmsecorrcoef30=mse(noisetarget30tr-noisea10tr)
valmsecorrcoef30=mse(noisetarget30tv-noisea10tv)
testmsecorrcoef30=mse(noisetarget30tst-noisea10tst)
 
trmsecorrcoef40=mse(noisetarget40tr-noisea10tr)
valmsecorrcoef40=mse(noisetarget40tv-noisea10tv)
testmsecorrcoef40=mse(noisetarget40tst-noisea10tst)
 
trmsecorrcoef50=mse(noisetarget50tr-noisea10tr)
valmsecorrcoef50=mse(noisetarget50tv-noisea10tv)
testmsecorrcoef50=mse(noisetarget50tst-noisea10tst)
 
r2msetr=[trmsecorrcoef10 trmsecorrcoef20 trmsecorrcoef30 trmsecorrcoef40 trmsecorrcoef50];
r2msetv=[valmsecorrcoef10 valmsecorrcoef20 valmsecorrcoef30 valmsecorrcoef40 valmsecorrcoef50];
r2msetst=[testmsecorrcoef10 testmsecorrcoef20 testmsecorrcoef30 testmsecorrcoef40 testmsecorrcoef50];
 
figure(71)
countr2mse=[r2msetr' r2msetv' r2msetst']
ybar2mse=countr2mse(1:5,:);
bar(ybar2mse,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('MSE Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2mse=mean(countr2mse')
%Input noise analysis based on R analysis
%%%%%%%%%% R analysis
in10tr=a10tr;
in20tr=a20tr;
in30tr=a30tr;
in40tr=a40tr;
in50tr=a50tr;
 
in10tv=a10v;
in20tv=a20v;
in30tv=a30v;
in40tv=a40v;
in50tv=a50v;
 
in10tst=a10t;
in20tst=a20t;
in30tst=a30t;
in40tst=a40t;
in50tst=a50t;
 
out10tr=tr;
out20tr=tr;
out30tr=tr;
out40tr=tr;
out50tr=tr;
 
out10tv=tav;
out20tv=tav;
out30tv=tav;
out40tv=tav;
out50tv=tav;
 
out10tst=tat;
out20tst=tat;
out30tst=tat;
out40tst=tat;
out50tst=tat;
 
trcorrcoef10=corrcoef(in10tr,out10tr)
valcorrcoef10=corrcoef(in10tv,out10tv)
testcorrcoef10=corrcoef(in10tst,out10tst)
 
noisetrcod10=diag(trcorrcoef10,1).*diag(trcorrcoef10,1)
noisevalcod10=diag(valcorrcoef10,1).*diag(valcorrcoef10,1)
noisetestcod10=diag(testcorrcoef10,1).*diag(testcorrcoef10,1)
%%%%
trcorrcoef20=corrcoef(in20tr,out20tr)
valcorrcoef20=corrcoef(in20tv,out20tv)
testcorrcoef20=corrcoef(in20tst,out20tst)
 
noisetrcod20=diag(trcorrcoef20,1).*diag(trcorrcoef20,1)
noisevalcod20=diag(valcorrcoef20,1).*diag(valcorrcoef20,1)
noisetestcod20=diag(testcorrcoef20,1).*diag(testcorrcoef20,1)
%%%%%%%%%%
trcorrcoef30=corrcoef(in30tr,out30tr)
valcorrcoef30=corrcoef(in30tv,out30tv)
testcorrcoef30=corrcoef(in30tst,out30tst)
 
noisetrcod30=diag(trcorrcoef30,1).*diag(trcorrcoef30,1)
noisevalcod30=diag(valcorrcoef30,1).*diag(valcorrcoef30,1)
noisetestcod30=diag(testcorrcoef30,1).*diag(testcorrcoef30,1)
%%%%%%%%%%%%%%%%%%%5
trcorrcoef40=corrcoef(in40tr,out40tr)
valcorrcoef40=corrcoef(in40tv,out40tv)
testcorrcoef40=corrcoef(in40tst,out40tst)
 
noisetrcod40=diag(trcorrcoef40,1).*diag(trcorrcoef40,1)
noisevalcod40=diag(valcorrcoef40,1).*diag(valcorrcoef40,1)
noisetestcod40=diag(testcorrcoef40,1).*diag(testcorrcoef40,1)
%%%%%%
trcorrcoef50=corrcoef(in50tr,out50tr)
valcorrcoef50=corrcoef(in50tv,out50tv)
testcorrcoef50=corrcoef(in50tst,out50tst)
 
noisetrcod50=diag(trcorrcoef50,1).*diag(trcorrcoef50,1)
noisevalcod50=diag(valcorrcoef50,1).*diag(valcorrcoef50,1);
noisetestcod50=diag(testcorrcoef50,1).*diag(testcorrcoef50,1)
 
r2tr=[noisetrcod10 noisetrcod20 noisetrcod30 noisetrcod40 noisetrcod50];
r2tv=[noisevalcod10 noisevalcod20 noisevalcod30 noisevalcod40 noisevalcod50];
r2tst=[noisetestcod10 noisetestcod20 noisetestcod30 noisetestcod40 noisetestcod50];
 
r2tr=sqrt(r2tr)
r2tv=sqrt(r2tv)
r2tst=sqrt(r2tst)
 
figure(72)
countr2=[r2tr' r2tv' r2tst']
ybar2=countr2(1:5,:);
bar(ybar2,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('R ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('R Analysis in Noisy Input Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2=mean(countr2')

.........................................................................................................................................
 
% % disp('Press any key to see a plot of the data together with the sine function.')
% % pause;
% 
clc;
clear all;
close all;
 
%Training data;
load gravdeep.dat;
aa=gravdeep(:,1);%longitude
bb=gravdeep(:,2);%latitude
cc=gravdeep(:,3);%BG value
dd=gravdeep(:,4);%Altitude
ee=gravdeep(:,5);%deep basement
% %----------------------------------------------------
%Training data;
% load gravshalow.dat;%205 no.
% a=gravshalow(:,1);%longitude
% b=gravshalow(:,2);%latitude
% c=gravshalow(:,3);%BG value
% d=gravshalow(:,4);%altitude
% e=gravshalow(:,5);%shallow basement
%%*************************TRAINING NOISY DATA****************************************
% aa=aa';
% bb=bb';
% cc=cc';
 
totalden = aa;
totalporo = bb;
totalgammar = cc;
%........................................
 
%Noise stability test
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
Aa1=0.99;%0.99
dsm(1)=a(1);
for i=2:nmax,
dsm(i)=Aa1.*dsm(i-1)+r(i);
end
dsm=dsm./max(abs(dsm));
dsm=dsm';
%r=r';
%den0=den.*r;
%den=a;
%totalden0=dsm;
%fp=0.0095;
totalden0=(totalden).*dsm;
totalden10=totalden+0.0010.*totalden0;
totalden20=totalden+0.0020.*totalden0;
totalden30=totalden+0.0030.*totalden0;
totalden40=totalden+0.0040.*totalden0;
totalden50=totalden+0.0050.*totalden0;
 
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
pr(1)=r(1);
Ba2=-0.0957;%-0.0957
for i=2:nmax,
pr(i)=Ba2.*pr(i-1)+r(i);
end
 
pr=pr./max(abs(pr));
pr=pr';
%poro0=poro.*r;
%poro=b;
%totalporo0=pr;
totalporo0=(totalporo).*pr;
totalporo10=totalporo+0.0010.*totalporo0;
totalporo20=totalporo+0.0020.*totalporo0;
totalporo30=totalporo+0.0030.*totalporo0;
totalporo40=totalporo+0.0040.*totalporo0;
totalporo50=totalporo+0.0050.*totalporo0;
 
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
gg(1)=r(1);
Ca3=0.0142;%0.0142
for i=2:nmax,
gg(i)=Ca3.*gg(i-1)+r(i);
end
gg=gg./max(abs(gg));
gg=gg';
%gammar0=gammar.*r;
%gammar=gammar';
%totalgammar0=gg;
totalgammar0=(totalgammar).*gg;
totalgammar10=totalgammar+0.0010.*totalgammar0;
totalgammar20=totalgammar+0.0020.*totalgammar0;
totalgammar30=totalgammar+0.0030.*totalgammar0;
totalgammar40=totalgammar+0.0040.*totalgammar0;
totalgammar50=totalgammar+0.0050.*totalgammar0;
 
fid = fopen('totalinerror10.dat','w');
totalerror10=[totalden10';totalporo10';totalgammar10'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror10);
fclose(fid);
 
fid = fopen('totalinerror20.dat','w');
totalerror20=[totalden20';totalporo20';totalgammar20'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror20);
fclose(fid);
 
fid = fopen('totalinerror30.dat','w');
totalerror30=[totalden30';totalporo30';totalgammar30'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror30);
fclose(fid);
 
fid = fopen('totalinerror40.dat','w');
totalerror40=[totalden40';totalporo40';totalgammar40'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror40);
fclose(fid);
 
fid = fopen('totalinerror50.dat','w');
totalerror50=[totalden50';totalporo50';totalgammar50'];
fprintf(fid,'%6.4f %6.4f %6.4f\n',totalerror50);
fclose(fid);
 
load totalinerror10.dat
load totalinerror20.dat
load totalinerror30.dat
load totalinerror40.dat
load totalinerror50.dat
 
%Training data;
totalinerror=totalinerror30;
aan=totalinerror(:,1);%longitude and n for nose
bbn=totalinerror(:,2);%latitude
ccn=totalinerror(:,3);%BG value
 
%........................................................................
target=ee;
 
randn('state',0);
nmax=211;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
Ta1=0.99;%0.99
dm(1)=a(1);
for i=2:nmax,
dm(i)=Ta1.*dm(i-1)+r(i);
end
 
dm=dm./max(abs(dm));
dm=dm';
%r=r';
%den0=den.*r;
%den=a;
%target0=dm';
target0=(target).*dm;
target10=target+0.10.*target0;
target20=target+0.20.*target0;
target30=target+0.30.*target0;
target40=target+0.40.*target0;
target50=target+0.50.*target0;
 
%************ANFIS*************************************************
% Real data
load gravwhole.dat;
ff=gravwhole(:,1);%longitude
gg=gravwhole(:,2);%latitude
hh=gravwhole(:,3);%BG value
ii=gravwhole(:,4);%altitude
%*********************************************************
 
x=[aan,bbn];
%t=ee;
t=target30;% target is corrupted with % white noise
x=x';
t=t';
 
%.....................................................noise adding......
 %ndata=211;
 noise = 0.15;          % Standard deviation of noise distribution.
 
 
%***********************************************************************
[xn,minp,maxp,tn,mint,maxt] = premnmx(x,t);
[xtrans,transMat] = prepca(xn,0.0000000000000002);
[R,Q] = size(xtrans);
iitst = 2:4:Q;
iival = 4:4:Q;
iitr = [1:4:Q 3:4:Q];
valP = xtrans(:,iival); val.T = tn(:,iival);
testP = xtrans(:,iitst); test.T = tn(:,iitst);
xtr = xtrans(:,iitr); ttr = tn(:,iitr);
 
xr=xtr';
zr=ttr';
 
trx=xr;
trt=zr;
 
trData=[trx,trt];
% x=[x1v,x2v];
 
% t=[x3v];
% xsm=[aa,bb,cc];
% tsm=ee;
 
%Data=[xsm,tsm];
 
xtrans=xtrans';
tn=tn';
 
Data=[xtrans,tn];
 
%numMFs=5;
%mfType='gbellmf';
%epoch_n=20;
%in_fis = genfis1(trnData,numMFs,mfType);
%out_fis = anfis(trnData,in_fis,20);
%evf=evalfis(x,out_fis);
%[afis] = postmnmx(evf,mint,maxt);
%plot(x1(1:350),x3(3:352),x1(1:350),afis(1:350));
%legend('Training Data','ANFIS Output');
%legend( 'Observed','ANFIS Network');
 
numMFs=2;
mfType='gbellmf';
 
trnData1=Data(1:4:Q, :);
trnData2=Data(3:4:Q, :);
 
valData=Data(4:4:Q, :);
testData=Data(2:4:Q, :);
 
trnData=[trnData1' trnData2'];
trnData=trnData';
 
 
%numMFs=5;
%mfType='gbellmf';
fismat = genfis1(trnData,numMFs,mfType);
%[fis1,error1,ss,fis2,error2] = anfis(trnData,options);
%fismat = genfis1(trnData);
 
%figure(15)
subplot(2,2,1)
plotmf(fismat, 'input', 1);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('Longitude','Fontsize',24);
ylabel('Degree of Membership','Fontsize',24);
subplot(2,2,2)
plotmf(fismat, 'input', 2)
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('Latitude','Fontsize',18);
ylabel('Degree of Membership','Fontsize',18);
subplot(2,2,3)
plotmf(fismat, 'input', 1)
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',20);axis tight;
xlabel('BGA','Fontsize',24);
ylabel('Degree of Membership','Fontsize',24);
subplot(2,2,4)
plotmf(fismat, 'input', 2)
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('BGA','Fontsize',18);
ylabel('Degree of Membership','Fontsize',24);
%[fismat1,error1,ss,fismat2,error2] = ...
   %anfis(trnData,fismat,[],[],chkData);
 
 
trnOpt(1)=500;
trnOpt(2)=0;
trnOpt(3)=.2;
trnOpt(4)=0.9;
trnOpt(5)=1.1;
dispOpt=1.0;
method=1;
[fismat1,trnError,ss,fismat2,valError] = ...
anfis(trnData,fismat,trnOpt,dispOpt,valData,method);
% [fismat1,trnError,ss,fismat2,valError,testError] = ...
% anfis(trnData,fismat,trnOpt,dispOpt,valData,testData,method);
 
%figure(16)
subplot(2,2,1)
plotmf(fismat2, 'input', 1)
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('Longitude','Fontsize',24);
ylabel('Degree of Membership','Fontsize',24);
subplot(2,2,2)
plotmf(fismat2, 'input', 2)
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('Latitude','Fontsize',18);
ylabel('Degree of Membership','Fontsize',24);
subplot(2,2,3)
plotmf(fismat2, 'input', 1)
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('Variable','Fontsize',24);
ylabel('Degree of Membership','Fontsize',24);
subplot(2,2,4)
plotmf(fismat2, 'input', 2)
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',24);axis tight;
xlabel('Lag distance','Fontsize',24);
ylabel('Degree of Membership','Fontsize',24);
 
%figure(17);
plot(trnError,'r','LineWidth',4);
xlabel('Epoch','Fontsize',24);
ylabel('MSE ','Fontsize',24);
title('Error Curves','Fontsize',24);
legend('Training Error');
grid on;
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
 
%figure(18);
plot(valError,'r','LineWidth',4);
xlabel('Epoch','Fontsize',24);
ylabel('MSE ','Fontsize',24);
title('Error Curves','Fontsize',24);
grid on;
legend('Testing Error');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',24);axis tight;
 
y1 = evalfis([trnData(:,1:2)],fismat2);
[y1] = postmnmx(y1',mint,maxt);
%anfis_output=y1;
% y2 = evalfis([chkData(:,1:2)],fismat2);
 
y2 = evalfis([valData(:,1:2)],fismat2);
[y2] = postmnmx(y2',mint,maxt);
 
y3 = evalfis([testData(:,1:2)],fismat2);
[y3] = postmnmx(y3',mint,maxt);
 
%ytro=[y1 y2];
 
 
xnew=x;
anfis_output = evalfis([xnew(:,1:2)],fismat2);
 
%Real data analysis
vinerror=[ff gg];
[e50v]= tramnmx(vinerror',minp,maxp);
[p2e50v] = trapca(e50v,transMat);
p2e50v=p2e50v';
arealv = evalfis([p2e50v(:,1:2)],fismat2);
%a50v = gpfwd(net,p2e50v');
[arealv] = postmnmx(arealv',mint,maxt);
 
% xreal=vinerror10;
%yreal = evalfis([xreal(:,1:3)],fismat2);
 
% trnT=t(1:4:Q, :);
% chkT=t(3:4:Q, :);
 
trnT1=tn(1:4:Q, :);
trnT2=tn(3:4:Q, :);
 
valT=tn(4:4:Q, :);
testT=tn(2:4:Q,:);
 
trnT=[trnT1' trnT2'];
trnT=trnT';
 
ee=ee';
tr=ee(:,iitr);
tav=ee(:,iival);
tat=ee(:,iitst);
 
%figure
 
subplot(3,1,1)
%plot(tr);hold on;plot(ypred)
plot(tr,'o-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6); hold on;plot(y1,'r-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
legend('Data','ANFISreg.');
%plot(tat);hold on;plot(ate);
%xlabel('No. of Data','Fontsize',18);
ylabel('Depth[Km]','Fontsize',24);
title('Training Interval:Original Vs. ANFISreg')
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
subplot(3,1,2)
%plot(tav);hold on;plot(av);
plot(tav,'o-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6); hold on;plot(y2,'r-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
legend('Data','ANFISreg.');
%plot(tat);hold on;plot(ate);
%xlabel('No. of Data','Fontsize',18);
ylabel('Depth[Km]','Fontsize',24);
title('Validation Interval:Original Vs. ANFISreg')
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
subplot(3,1,3)
plot(tat,'o-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6); hold on;plot(y3,'r-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
legend('Data','ANFISreg.');
%plot(tat);hold on;plot(ate);
xlabel('No. of Data','Fontsize',24);
ylabel('Depth[Km]','Fontsize',24);
title('Test Interval:Oroginal Vs. ANFISreg');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
%figure
plotregression(tr,y1,'Regression:Training Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
xlabel('Target Depth[Km]','Fontsize',24);
ylabel('Predicted Depth[Km]','Fontsize',24);
 
%legend('Training Interval')
dbp=tr-y1;
msedbp=mse(dbp)%mse of trainind set
maedbp=mae(dbp)
retrain=1.0-sum((tr-y1).^2)./sum(tr.^2)
dtrain=1.0-sum((tr-y1).^2)./sum((y1-mean(tr)+tr-mean(tr)).^2)
Rtrain=corrcoef(tr,y1)
 
% figure
%subplot(3,2,4)
plotregression(tav,y2,'Regression:Validation Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
xlabel('Target Depth[Km]','Fontsize',24);
ylabel('Predicted Depth[Km]','Fontsize',24);
 
%legend('Validation Interval')
vdbp=tav-y2;
msevdbp=mse(vdbp)%mse of validation set
maevdbp=mae(vdbp)
reval=1.0-sum((tav-y2).^2)./sum(tav.^2)
dval=1.0-sum((tav-y2).^2)./sum((y2-mean(y2)+tav-mean(tav)).^2)
Rval=corrcoef(tav,y2)
% figure
%subplot(3,2,6)
plotregression(tat,y3,'Regression:Test Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
xlabel('Target Depth[Km]','Fontsize',24);
ylabel('Predicted Depth[Km]','Fontsize',24);
 
tdbp=tat-y3;
msetdbp=mse(tdbp)%mse of test set
maetdbp=mae(tdbp)
retest=1.0-sum((tat-y3).^2)./sum(tat.^2)
dtest=1.0-sum((tat-y3).^2)./sum((y3-mean(y3)+tat-mean(tat)).^2)
Rtest=corrcoef(tat,y3)
 
trxis=1:106;
[trp,trs]=polyfit(trxis,y1,1);
[trp_sdepthtr,deltatr]=polyval(trp,trxis,trs);
trsig=2.0.*deltatr;
 
valxis=1:52;
[valp,vals]=polyfit(valxis,y2,1);
[valp_sdepthval,deltaval]=polyval(valp,valxis,vals);
valsig=2.0.*deltaval;
 
tstxis=1:53;
[tstp,tsts]=polyfit(tstxis,y3,1);
[tstp_sdepthtst,deltatst]=polyval(tstp,tstxis,tsts);
tstsig=2.0.*deltatst;
%..........................................................................
% 
% 
 
% 
 
 
figure(26);
 
xis=1:106;
E1=trsig.*ones(size(y1));
% %subplot(3,1,1),
errorbar(xis,y1,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tr,'ro-',xis, y1,'o-','LineWidth',4,...
                 'MarkerEdgeColor','k',...
                 'MarkerFaceColor','g',...
                 'MarkerSize',6);grid on;axis([0 106 -6 6]);
 legend('Error-bar','Target','ANFISreg.');
 xlabel('No. of Data','Fontsize',24);
 ylabel('Sediment Depth(Km)','Fontsize',24);
 title('Training Interval')
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
figure(27)
 
xis=1:52;
E1=valsig.*ones(size(y2));
% %subplot(3,1,1),
errorbar(xis,y2,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tav,'ro-',xis, y2,'o-','LineWidth',4,...
                 'MarkerEdgeColor','k',...
                 'MarkerFaceColor','g',...
                 'MarkerSize',6);grid on;axis([0 106 -6 6]);
 legend('Error-bar','Target','ANFISreg.');
 xlabel('No. of Data','Fontsize',24);
 ylabel('Sediment Depth(Km)','Fontsize',24);
 title('Validation Interval')
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
figure(28)
 
xis=1:53;
E1=tstsig.*ones(size(y3));
% %subplot(3,1,1),
errorbar(xis,y3,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tat,'ro-',xis, y3,'o-','LineWidth',4,...
                 'MarkerEdgeColor','k',...
                 'MarkerFaceColor','g',...
                 'MarkerSize',6);grid on;axis([0 106 -6 6]);
 legend('Error-bar','Target','ANFISreg.');
 xlabel('No. of Data','Fontsize',24);
 ylabel('Sediment Depth(Km)','Fontsize',24);
 title('Test Interval')
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%........................................
 
%clear all
%Statistical %Training....................................................
%.........................................................................
ntotal=211;
ntrain=106;
nval=52;
ntest=53;
 
tr=ee(:,iitr);
tav=ee(:,iival);
tat=ee(:,iitst);
 
x1train=tr(1:ntrain);
x1val=tav(1:nval);
x1test=tat(1:ntest);
 
x2train=y1(1:ntrain);
x2val=y2(1:nval);
x2test=y3(1:ntest);
 
disp('Statistics of x1/obs......................................');
x1avgtrain=mean(x1train)
x1avgnval=mean(x1val)
x1avgtest=mean(x1test)
 
x1stdtrain=std(x1train)
x1stdval=std(x1val)
x1stdtest=std(x1test)
 
x1trmedian=median(x1train)
x1valmedian=median(x1val)
x1testmedian=median(x1test)
 
x1trmode=mode(x1train)
x1valmode=mode(x1val)
x1testmode=mode(x1test)
 
x1trskewness=skewness(x1train)
x1valskewness=skewness(x1val)
x1testskewness=skewness(x1test)
 
x1trkurtosis=kurtosis(x1train)
x1valkurtosis=kurtosis(x1val)
x1testkurtosis=kurtosis(x1test)
 
disp('Statistics of x2/ANFISreg......................................')
x2avgtrain=mean(x2train)
x2avgnval=mean(x2val)
x2avgtest=mean(x2test)
 
x2stdtrain=std(x2train)
x2stdval=std(x2val)
x2stdtest=std(x2test)
 
x2trmedian=median(x2train);
x2valmedian=median(x2val)
x2testmedian=median(x2test)
 
x2trmode=mode(x2train)
x2valmode=mode(x2val)
x2testmode=mode(x2test)
 
x2trskewness=skewness(x2train)
x2valskewness=skewness(x2val)
x2testskewness=skewness(x2test)
 
x2trkurtosis=kurtosis(x2train)
x2valkurtosis=kurtosis(x2val)
x2testkurtosis=kurtosis(x2test)
 
x2errtrain=mse(x2train-x1train)
x2errval=mse(x2val-x1val)
x2errtest=mse(x2test-x1test)
 
x2retrain=1.0-sum((x1train-x2train).^2)./sum(x1train.^2)
x2reval=1.0-sum((x1val-x2val).^2)./sum(x1val.^2)
x2retest=1.0-sum((x1test-x2test).^2)./sum(x1test.^2)
 
x2dtrain=1.0-sum((x1train-x2train).^2)./sum(abs(x2train-mean(x1train))+abs(x1train-mean(x1train)).^2)
x2dval=1.0-sum((x1val-x2val).^2)./sum(abs(x2val-mean(x1val))+abs(x1val-mean(x1val)).^2)
x2dtest=1.0-sum((x1test-x2test).^2)./sum(abs(x2test-mean(x1test))+abs(x1test-mean(x1test)).^2)
 
trcorrcoef=corrcoef(x1train,x2train)
valcorrcoef=corrcoef(x1val,x2val)
testcorrcoef=corrcoef(x1test,x2test)
 
trcod=diag(trcorrcoef,1).*diag(trcorrcoef,1)
valcod=diag(valcorrcoef,1).*diag(valcorrcoef,1)
testcod=diag(testcorrcoef,1).*diag(testcorrcoef,1)
 
%***************Training Error Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error analysis**********************************************************************
aa=aa';
bb=bb';
cc=cc';
% inputden = xr(:,1);
% inputporo = xr(:,2);
% inputgammar = xr(:,3);
 
inputden = aa(:,iitr);
inputporo = bb(:,iitr);
inputgammar = cc(:,iitr);
 
inputden=inputden';
inputporo=inputporo';
inputgammar=inputgammar';
 
randn('state',0);
nmax=106;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
a(1)=r(1);
a1=0.94;%0.94
dtr(1)=a(1);
for i=2:nmax,
dtr(i)=a1.*dtr(i-1)+r(i);
end
 
dtr=dtr./max(abs(dtr));
dtr=dtr';
%r=r';
%den0=den.*r;
%den=a;
%%inputden0=dtr;
inputden0=(inputden).*dtr;
inputden10=inputden+0.0010.*inputden0;
inputden20=inputden+0.0020.*inputden0;
inputden30=inputden+0.0030.*inputden0;
inputden40=inputden+0.0040.*inputden0;
inputden50=inputden+0.0050.*inputden0;
 
randn('state',0);
nmax=106;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
ptr(1)=r(1);
a2=-0.107;%-0.107
for i=2:nmax,
ptr(i)=a2.*ptr(i-1)+r(i);
end
ptr=ptr';
ptr=ptr./max(abs(ptr));
%poro0=poro.*r;
%poro=b;
%%inputporo0=ptr;
inputporo0=(inputporo).*ptr;
inputporo10=inputporo+0.0010.*inputporo0;
inputporo20=inputporo+0.0020.*inputporo0;
inputporo30=inputporo+0.0030.*inputporo0;
inputporo40=inputporo+0.0040.*inputporo0;
inputporo50=inputporo+0.0050.*inputporo0;
 
randn('state',0);
nmax=106;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
ggtr(1)=r(1);
a3=0.122;%0.122
for i=2:nmax,
ggtr(i)=a3.*ggtr(i-1)+r(i);
end
 
ggtr=ggtr./max(abs(ggtr));
ggtr=ggtr';
%gammar0=gammar.*r;
%gammar=gammar';
%%inputgammar0=ggtr;
inputgammar0=(inputgammar).*ggtr;
inputgammar10=inputgammar+0.0010.*inputgammar0;
inputgammar20=inputgammar+0.0020.*inputgammar0;
inputgammar30=inputgammar+0.0030.*inputgammar0;
inputgammar40=inputgammar+0.0040.*inputgammar0;
inputgammar50=inputgammar+0.0050.*inputgammar0;
 
fid = fopen('sinerror10.dat','w');
%inputerror10=[inputden10';inputporo10';inputgammar10'];
inputerror10=[inputden10';inputporo10'];
fprintf(fid,'%6.4f %6.4f\n',inputerror10);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror10);
fclose(fid);
 
fid = fopen('sinerror20.dat','w');
%inputerror20=[inputden20';inputporo20';inputgammar20'];
inputerror20=[inputden20';inputporo20'];
fprintf(fid,'%6.4f %6.4f\n',inputerror20);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror20);
fclose(fid);
 
fid = fopen('sinerror30.dat','w');
%inputerror30=[inputden30';inputporo30';inputgammar30'];
inputerror30=[inputden30';inputporo30'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror30);
fprintf(fid,'%6.4f %6.4f\n',inputerror30);
fclose(fid);
 
fid = fopen('sinerror40.dat','w');
%inputerror40=[inputden40';inputporo40';inputgammar40'];
inputerror40=[inputden40';inputporo40'];
fprintf(fid,'%6.4f %6.4f\n',inputerror40);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror40);
fclose(fid);
 
fid = fopen('sinerror50.dat','w');
%inputerror50=[inputden50';inputporo50';inputgammar50'];
inputerror50=[inputden50';inputporo50'];
fprintf(fid,'%6.4f %6.4f\n',inputerror50);
%fprintf(fid,'%6.4f %6.4f %6.4f\n',inputerror50);
fclose(fid);
 
load sinerror10.dat;
load sinerror20.dat;
load sinerror30.dat;
load sinerror40.dat;
load sinerror50.dat
 
 
[e10tr]= tramnmx(sinerror10',minp,maxp);
[p2e10tr] = trapca(e10tr,transMat);
p2e10tr=p2e10tr';
a10tr = evalfis([p2e10tr(:,1:2)],fismat2);
%a10r = gpfwd(net,p2e10');
[a10tr] = postmnmx(a10tr',mint,maxt);
 
 
 
%t=t(:,iitr);
%tr=t(:,iitr);
tr=tr';
trdev10=tr-a10tr;
trdev10=trdev10';
trdp10=trdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
trmse10=mse(trdp10);
 
figure(29);
 
subplot(1,1,1),plot(trdp10,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e20tr]= tramnmx(sinerror20',minp,maxp);
[p2e20tr] = trapca(e20tr,transMat);
p2e20tr=p2e20tr';
a20tr = evalfis([p2e20tr(:,1:2)],fismat2);
%a20tr = gpfwd(net,p2e20tr');
[a20tr] = postmnmx(a20tr',mint,maxt);
 
 
 
trdev20=tr-a20tr;
trdev20=trdev20';
trdp20=trdev20(:,1);
trmse20=mse(trdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(30);
 
subplot(1,1,1),plot(trdp20,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e30tr]= tramnmx(sinerror30',minp,maxp);
[p2e30tr] = trapca(e30tr,transMat);
p2e30tr=p2e30tr';
a30tr = evalfis([p2e30tr(:,1:2)],fismat2);
%a30tr = gpfwd(net,p2e30tr');
[a30tr] = postmnmx(a30tr',mint,maxt);
 
 
 
trdev30=tr-a30tr;
trdev30=trdev30';
trdp30=trdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
trmse30=mse(trdp30);
 
figure(31);
 
subplot(1,1,1),plot(trdp30,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e40tr]= tramnmx(sinerror40',minp,maxp);
[p2e40tr] = trapca(e40tr,transMat);
p2e40tr=p2e40tr';
a40tr = evalfis([p2e40tr(:,1:2)],fismat2);
%a40tr = gpfwd(net,p2e40tr');
[a40tr] = postmnmx(a40tr',mint,maxt);
 
trdev40=tr-a40tr;
trdev40=trdev40';
trdp40=trdev40(:,1);
trmse40=mse(trdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(32);
 
subplot(1,1,1),plot(trdp40,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e50tr]= tramnmx(sinerror50',minp,maxp);
[p2e50tr] = trapca(e50tr,transMat);
p2e50tr=p2e50tr';
a50tr = evalfis([p2e50tr(:,1:2)],fismat2);
%a50tr = gpfwd(net,p2e50tr');
[a50tr] = postmnmx(a50tr',mint,maxt);
 
 
 
trdev50=tr-a50tr;
trdev50=trdev50';
trdp50=trdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
trmse50=mse(trdp50);
 
figure(33);
subplot(1,1,1),plot(trdp50,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Training Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
%error diagram
%figure(34)
 
nx=[10 20 30 40 50];
msetr=[trmse10 trmse20 trmse30 trmse40 trmse50]
bar(nx,msetr);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Training Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Obs.','ANFISreg', 'ANFISreg+STD','ANFISreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%******************Validation Error ANALYSIS****************************
 
% vdenv = valP(:,1);
% vporov = valP(:,2);
% vgammarv = valP(:,3);
 
vdenv = aa(:,iival);
vporov = bb(:,iival);
vgammarv = cc(:,iival);
 
valden=vdenv';
valporo=vporov';
valgammar=vgammarv';
 
randn('state',0);
nmax=52;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
dv(1)=r(1);
A1=0.99;%0.99
for i=2:nmax,
dv(i)=A1.*dv(i-1)+r(i);
end
 
dv=dv./max(abs(dv));
dv=dv';
%r=r';
%den0=den.*r;
%%valden0=dv;
valden0=(valden).*dv;
valden10=valden+0.0010.*valden0;
valden20=valden+0.0020.*valden0;
valden30=valden+0.0030.*valden0;
valden40=valden+0.0040.*valden0;
valden50=valden+0.0050.*valden0;
 
randn('state',0);
nmax=52;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
pv(1)=r(1);
A2=-0.0495;%-0.0495
for i=2:nmax,
pv(i)=A2.*pv(i-1)+r(i);
end
pv=pv';
pv=pv./max(abs(pv));
%poro0=poro.*r;
%poro=b;
%valporo0=pv;
valporo0=(valporo).*pv;
valporo10=valporo+0.0010.*valporo0;
valporo20=valporo+0.0020.*valporo0;
valporo30=valporo+0.0030.*valporo0;
valporo40=valporo+0.0040.*valporo0;
valporo50=valporo+0.0050.*valporo0;
 
randn('state',0);
nmax=52;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
gva(1)=r(1);
A3=0.1408;%0.1408
for i=2:nmax,
gva(i)=A3.*gva(i-1)+r(i);
end
 
gva=gva./max(abs(gva));
gva=gva';
%gammar0=gammar.*r;
%gammar=gammar';
%valgammar0=gva;
valgammar0=(valgammar).*gva;
valgammar10=valgammar+0.0010.*valgammar0;
valgammar20=valgammar+0.0020.*valgammar0;
valgammar30=valgammar+0.0030.*valgammar0;
valgammar40=valgammar+0.0040.*valgammar0;
valgammar50=valgammar+0.0050.*valgammar0;
 
fid = fopen('vinerror10.dat','w');
%valerror10=[valden10';valporo10';valgammar10'];
valerror10=[valden10';valporo10'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror10);
fprintf(fid,'%6.4f %6.4f\n',valerror10);
fclose(fid);
 
fid = fopen('vinerror20.dat','w');
%valerror20=[valden20';valporo20';valgammar20'];
valerror20=[valden20';valporo20'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror20);
fprintf(fid,'%6.4f %6.4f\n',valerror20);
fclose(fid);
 
fid = fopen('vinerror30.dat','w');
%valerror30=[valden30';valporo30';valgammar30'];
valerror30=[valden30';valporo30'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror30);
fprintf(fid,'%6.4f %6.4f\n',valerror30);
fclose(fid);
 
fid = fopen('vinerror40.dat','w');
%valerror40=[valden40';valporo40';valgammar40'];
valerror40=[valden40';valporo40'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror40);
fprintf(fid,'%6.4f %6.4f\n',valerror40);
fclose(fid);
 
fid = fopen('vinerror50.dat','w');
%valerror50=[valden50';valporo50';valgammar50'];
valerror50=[valden50';valporo50'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',valerror50);
fprintf(fid,'%6.4f %6.4f\n',valerror50);
fclose(fid);
 
load vinerror10.dat;
load vinerror20.dat;
load vinerror30.dat;
load vinerror40.dat;
load vinerror50.dat
 
 
[e10v]= tramnmx(vinerror10',minp,maxp);
[p2e10v] = trapca(e10v,transMat);
p2e10v=p2e10v';
a10v = evalfis([p2e10v(:,1:2)],fismat2);
%a10v = gpfwd(net,p2e10v');
[a10v] = postmnmx(a10v',mint,maxt);
 
%tav=tav';
vdev10=tav'-a10v;
vdev10=vdev10';
vdp10=vdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
valmse10=mse(vdp10);
 
figure(35);
 
subplot(1,1,1),plot(vdp10,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e20v]= tramnmx(vinerror20',minp,maxp);
[p2e20v] = trapca(e20v,transMat);
p2e20v=p2e20v';
a20v = evalfis([p2e20v(:,1:2)],fismat2);
%a20v = gpfwd(net,p2e20v');
[a20v] = postmnmx(a20v',mint,maxt);
 
vdev20=tav'-a20v;
vdev20=vdev20';
vdp20=vdev20(:,1);
valmse20=mse(vdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(36);
 
subplot(1,1,1),plot(vdp20,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e30v]= tramnmx(vinerror30',minp,maxp);
[p2e30v] = trapca(e30v,transMat);
p2e30v=p2e30v';
a30v = evalfis([p2e30v(:,1:2)],fismat2);
%a30v = gpfwd(net,p2e30v');
[a30v] = postmnmx(a30v',mint,maxt);
 
 
 
vdev30=tav'-a30v;
vdev30=vdev30';
vdp30=vdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
valmse30=mse(vdp30);
 
figure(37);
 
subplot(1,1,1),plot(vdp30,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e40v]= tramnmx(vinerror40',minp,maxp);
[p2e40v] = trapca(e40v,transMat);
p2e40v=p2e40v';
a40v = evalfis([p2e40v(:,1:2)],fismat2);
%a40v = gpfwd(net,p2e40v');
[a40v] = postmnmx(a40v',mint,maxt);
 
vdev40=tav'-a40v;
vdev40=vdev40';
vdp40=vdev40(:,1);
valmse40=mse(vdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(38);
 
subplot(1,1,1),plot(vdp40,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e50v]= tramnmx(vinerror50',minp,maxp);
[p2e50v] = trapca(e50v,transMat);
p2e50v=p2e50v';
a50v = evalfis([p2e50v(:,1:2)],fismat2);
%a50v = gpfwd(net,p2e50v');
[a50v] = postmnmx(a50v',mint,maxt);
 
 
 
vdev50=tav'-a50v;
vdev50=vdev50';
vdp50=vdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
valmse50=mse(vdp50);
 
figure(40);
subplot(1,1,1),plot(vdp50,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Validation Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',24);axis tight
 
%error diagram
figure(41)
nx=[10 20 30 40 50];
mseval=[valmse10 valmse20 valmse30 valmse40 valmse50];
bar(nx,mseval);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Validation Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Obs.','ANFISreg', 'ANFISreg+STD','ANFISreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%**********************Test Error Analysis*****************************
%Error Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error analysis
 
% dent = testP(:,1);
% porot = testP(:,2);
% gammart = testP(:,3);
 
dent = aa(:,iitst);
porot = bb(:,iitst);
gammart = cc(:,iitst);
 
testden=dent';
testporo=porot';
testgammar=gammart';
 
randn('state',0);
nmax=53;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
dtst(1)=r(1);
B1=0.99;%0.99
 
for i=2:nmax,
dtst(i)=B1.*dtst(i-1)+r(i);
end
 
dtst=dtst./max(abs(dtst));
dtst=dtst';
%r=r';
%den0=den.*r;
%den=a;
%testden0=dtst;
testden0=(testden).*dtst;
testden10=testden+0.0010.*testden0;
testden20=testden+0.0020.*testden0;
testden30=testden+0.0030.*testden0;
testden40=testden+0.0040.*testden0;
testden50=testden+0.0050.*testden0;
 
randn('state',0);
nmax=53;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
pt(1)=r(1);
B2=-0.2317;%-0.2317
for i=2:nmax,
pt(i)=B2.*pt(i-1)+r(i);
end
 
pt=pt./max(abs(pt));
pt=pt';
%poro0=poro.*r;
%poro=b;
%testporo0=pt;
testporo0=(testporo).*pt;
testporo10=testporo+0.0010.*testporo0;
testporo20=testporo+0.0020.*testporo0;
testporo30=testporo+0.0030.*testporo0;
testporo40=testporo+0.0040.*testporo0;
testporo50=testporo+0.0050.*testporo0;
 
randn('state',0);
nmax=53;
% White time series
r=rand(1,nmax);
r=r/max(abs(r));
gt(1)=r(1);
B3=-0.036317;%-0.036317
for i=2:nmax,
gt(i)=B3.*gt(i-1)+r(i);
end
 
gt=gt./max(abs(gt));
gt=gt';
%gammar0=gammar.*r;
%gammar=gammar';
%testgammar0=gt;
testgammar0=(testgammar).*gt;
testgammar10=testgammar+0.0010.*testgammar0;
testgammar20=testgammar+0.0020.*testgammar0;
testgammar30=testgammar+0.0030.*testgammar0;
testgammar40=testgammar+0.0040.*testgammar0;
testgammar50=testgammar+0.0050.*testgammar0;
 
fid = fopen('tinerror10.dat','w');
%testerror10=[testden10';testporo10';testgammar10'];
testerror10=[testden10';testporo10'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror10);
fprintf(fid,'%6.4f %6.4f\n',testerror10);
fclose(fid);
 
fid = fopen('tinerror20.dat','w');
%testerror20=[testden20';testporo20';testgammar20'];
testerror20=[testden20';testporo20'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror20);
fprintf(fid,'%6.4f %6.4f\n',testerror20);
fclose(fid);
 
fid = fopen('tinerror30.dat','w');
%testerror30=[testden30';testporo30';testgammar30'];
testerror30=[testden30';testporo30'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror30);
fprintf(fid,'%6.4f %6.4f\n',testerror30);
fclose(fid);
 
fid = fopen('tinerror40.dat','w');
%testerror40=[testden40';testporo40';testgammar40'];
testerror40=[testden40';testporo40'];
%fprintf(fid,'%6.4f %6.4f %6.4f\n',testerror40);
fprintf(fid,'%6.4f %6.4f\n',testerror40);
fclose(fid);
 
fid = fopen('tinerror50.dat','w');
%testerror50=[testden50';testporo50';testgammar50'];
testerror50=[testden50';testporo50'];
%fprintf(fid,'%6.4f 6.4f %6.4f\n',testerror50);
fprintf(fid,'%6.4f %6.4f\n',testerror50);
fclose(fid);
 
load tinerror10.dat;
load tinerror20.dat;
load tinerror30.dat;
load tinerror40.dat;
load tinerror50.dat
 
 
[e10t]= tramnmx(tinerror10',minp,maxp);
[p2e10t] = trapca(e10t,transMat);
p2e10t=p2e10t';
a10t = evalfis([p2e10t(:,1:2)],fismat2);
%a10t = gpfwd(net,p2e10t');
[a10t] = postmnmx(a10t',mint,maxt);
 
 
%tat=tat';
tdev10=tat'-a10t;
tdev10=tdev10';
tdp10=tdev10(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
tstmse10=mse(tdp10);
 
figure(45);
 
subplot(1,1,1),plot(tdp10,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 10% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
[e20t]= tramnmx(tinerror20',minp,maxp);
[p2e20t] = trapca(e20t,transMat);
p2e20t=p2e20t';
a20t = evalfis([p2e20t(:,1:2)],fismat2);
%a20t = gpfwd(net,p2e20t');
[a20t] = postmnmx(a20t',mint,maxt);
 
 
 
tdev20=tat'-a20t;
tdev20=tdev20';
tdp20=tdev20(:,1);
tstmse20=mse(tdp20);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(46);
 
subplot(1,1,1),plot(tdp20,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 20% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
[e30t]= tramnmx(tinerror30',minp,maxp);
[p2e30t] = trapca(e30t,transMat);
p2e30t=p2e30t';
a30t = evalfis([p2e30t(:,1:2)],fismat2);
%a30t = gpfwd(net,p2e10t');
[a30t] = postmnmx(a30t',mint,maxt);
 
 
 
tdev30=tat'-a30t;
tdev30=tdev30';
tdp30=tdev30(:,1);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
tstmse30=mse(tdp30);
 
figure(47);
 
subplot(1,1,1),plot(tdp30,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 30% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
[e40t]= tramnmx(tinerror40',minp,maxp);
[p2e40t] = trapca(e40t,transMat);
p2e40t=p2e40t';
a40t = evalfis([p2e40t(:,1:2)],fismat2);
%a40t = gpfwd(net,p2e40t');
[a40t] = postmnmx(a40t',mint,maxt);
 
tdev40=tat'-a40t;
tdev40=tdev40';
tdp40=tdev40(:,1);
tstmse40=mse(tdp40);
% dm10=dev10(:,2);
% dh10=dev10(:,3);
 
figure(48);
 
subplot(1,1,1),plot(tdp40,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 40% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',24);axis tight
 
 
[e50t]= tramnmx(tinerror50',minp,maxp);
[p2e50t] = trapca(e50t,transMat);
p2e50t=p2e50t';
a50t = evalfis([p2e50t(:,1:2)],fismat2);
%a50t = gpfwd(net,p2e50t');
[a50t] = postmnmx(a50t',mint,maxt);
 
 
 
tdev50=tat'-a50t;
tdev50=tdev50';
tdp50=tdev50(:,1);
%dm10=dev10(:,2);
%dh10=dev10(:,3);
tstmse50=mse(tdp50);
 
figure(50);
 
subplot(1,1,1),plot(tdp50,'ro-','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -1 1]);
title('Error Deviation When Test Set Data is Corrupted with 50% Correlated Noise')
xlabel('No. of Data','Fontsize',24);
ylabel('Error Dev.','Fontsize',24);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
 
 
 
%figure(51)
nx=[10 20 30 40 50];
msetst=[tstmse10 tstmse20 tstmse30 tstmse40 tstmse50];
bar(nx,msetst);
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('Test Interval');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Obs.','ANFISreg', 'ANFISreg+STD','ANFISreg-STD');
%set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
figure(52)
count=[msetr' mseval' msetst']
ybar=count(1:5,:);
bar(ybar,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('MSE Analysis in Noisy Input Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
meancount=mean(count')
 
figure(53)
count=[mseval' msetst'];
ybar=count(1:5,:);
bar(ybar,'LineWidth',3,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',20);ylabel('MSE ','FontName','Tahoma','Fontsize',20);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',20);
% set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
title('MSE Analysis in Noisy Input Data');
legend('Validation','Test','Location','northeast');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
%..................................Regression 
% 
 
%Validation error...................................................
figure(60)
plotregression(tav,a10v)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
figure(61)
plotregression(tav,a20v)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
figure(62)
plotregression(tav,a30v)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
figure(63)
plotregression(tav,a40v)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
figure(64)
plotregression(tav,a50v)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Validation Data+10% Noise','Location','northeast');
figure(65)
plotregression(tat,a10t)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
figure(66)
plotregression(tat,a20t)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
figure(67)
plotregression(tat,a30t)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
figure(68)
plotregression(tat,a40t)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
figure(69)
plotregression(tat,a50t)
xlabel('Target','FontName','Tahoma','Fontsize',24);ylabel('ANFISreg','FontName','Tahoma','Fontsize',24);
%title('Input+10%Noise','FontName','Tahoma','Fontsize',20);
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
%legend('Test Data+10% Noise','Location','northeast');
 
%correlation analysis at different target noise
vinerror10=[ff,gg];
[e10]= tramnmx(vinerror10',minp,maxp);
[p2e10] = trapca(e10,transMat);
p2e10=p2e10';
a10 = evalfis([p2e10(:,1:2)],fismat2);
%a10v = gpfwd(net,p2e10v');
[a10] = postmnmx(a10',mint,maxt);
% a10 = gpfwd(net,p2e10');
% [a10] = postmnmx(a10',mint,maxt);
 
target10=target10';
target20=target20';
target30=target30';
target40=target40';
target50=target50';
 
noisetarget10tr=target10(:,iitr);
noisetarget20tr=target20(:,iitr);
noisetarget30tr=target30(:,iitr);
noisetarget40tr=target40(:,iitr);
noisetarget50tr=target50(:,iitr);
 
noisetarget10tv=target10(:,iival);
noisetarget20tv=target20(:,iival);
noisetarget30tv=target30(:,iival);
noisetarget40tv=target40(:,iival);
noisetarget50tv=target50(:,iival);
 
noisetarget10tst=target10(:,iitst);
noisetarget20tst=target20(:,iitst);
noisetarget30tst=target30(:,iitst);
noisetarget40tst=target40(:,iitst);
noisetarget50tst=target50(:,iitst);
 
noisea10tr=a10(:,iitr);
noisea10tv=a10(:,iival);
noisea10tst=a10(:,iitst);
 
%%%%%%%%%% R analysis
trcorrcoef10=corrcoef(noisetarget10tr,noisea10tr)
valcorrcoef10=corrcoef(noisetarget10tv,noisea10tv)
testcorrcoef10=corrcoef(noisetarget10tst,noisea10tst)
 
noisetrcod10=diag(trcorrcoef10,1).*diag(trcorrcoef10,1)
noisevalcod10=diag(valcorrcoef10,1).*diag(valcorrcoef10,1)
noisetestcod10=diag(testcorrcoef10,1).*diag(testcorrcoef10,1)
%%%%
trcorrcoef20=corrcoef(noisetarget20tr,noisea10tr)
valcorrcoef20=corrcoef(noisetarget20tv,noisea10tv)
testcorrcoef20=corrcoef(noisetarget20tst,noisea10tst)
 
noisetrcod20=diag(trcorrcoef20,1).*diag(trcorrcoef20,1)
noisevalcod20=diag(valcorrcoef20,1).*diag(valcorrcoef20,1)
noisetestcod20=diag(testcorrcoef20,1).*diag(testcorrcoef20,1)
%%%%%%%%%%
trcorrcoef30=corrcoef(noisetarget30tr,noisea10tr)
valcorrcoef30=corrcoef(noisetarget30tv,noisea10tv)
testcorrcoef30=corrcoef(noisetarget30tst,noisea10tst)
 
noisetrcod30=diag(trcorrcoef30,1).*diag(trcorrcoef30,1)
noisevalcod30=diag(valcorrcoef30,1).*diag(valcorrcoef30,1)
noisetestcod30=diag(testcorrcoef30,1).*diag(testcorrcoef30,1)
%%%%%%%%%%%%%%%%%%%5
trcorrcoef40=corrcoef(noisetarget40tr,noisea10tr)
valcorrcoef40=corrcoef(noisetarget40tv,noisea10tv)
testcorrcoef40=corrcoef(noisetarget40tst,noisea10tst)
 
noisetrcod40=diag(trcorrcoef40,1).*diag(trcorrcoef40,1)
noisevalcod40=diag(valcorrcoef40,1).*diag(valcorrcoef40,1)
noisetestcod40=diag(testcorrcoef40,1).*diag(testcorrcoef40,1)
%%%%%%
trcorrcoef50=corrcoef(noisetarget50tr,noisea10tr)
valcorrcoef50=corrcoef(noisetarget50tv,noisea10tv)
testcorrcoef50=corrcoef(noisetarget50tst,noisea10tst)
 
noisetrcod50=diag(trcorrcoef50,1).*diag(trcorrcoef50,1)
noisevalcod50=diag(valcorrcoef50,1).*diag(valcorrcoef50,1)
noisetestcod50=diag(testcorrcoef50,1).*diag(testcorrcoef50,1)
 
r2tr=[noisetrcod10 noisetrcod20 noisetrcod30 noisetrcod40 noisetrcod50];
r2tv=[noisevalcod10 noisevalcod20 noisevalcod30 noisevalcod40 noisevalcod50];
r2tst=[noisetestcod10 noisetestcod20 noisetestcod30 noisetestcod40 noisetestcod50];
 
r2tr=sqrt(r2tr)
r2tv=sqrt(r2tv)
r2tst=sqrt(r2tst)
 
figure(70)
countr2=[r2tr' r2tv' r2tst']
ybar2=countr2(1:5,:);
bar(ybar2,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('R ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('R Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2=mean(countr2')
%..................MSE analysis
trmsecorrcoef10=mse(noisetarget10tr-noisea10tr)
valmsecorrcoef10=mse(noisetarget10tv-noisea10tv)
testmsecorrcoef10=mse(noisetarget10tst-noisea10tst)
 
trmsecorrcoef20=mse(noisetarget20tr-noisea10tr)
valmsecorrcoef20=mse(noisetarget20tv-noisea10tv)
testmsecorrcoef20=mse(noisetarget20tst-noisea10tst)
 
trmsecorrcoef30=mse(noisetarget30tr-noisea10tr)
valmsecorrcoef30=mse(noisetarget30tv-noisea10tv)
testmsecorrcoef30=mse(noisetarget30tst-noisea10tst)
 
trmsecorrcoef40=mse(noisetarget40tr-noisea10tr)
valmsecorrcoef40=mse(noisetarget40tv-noisea10tv)
testmsecorrcoef40=mse(noisetarget40tst-noisea10tst)
 
trmsecorrcoef50=mse(noisetarget50tr-noisea10tr)
valmsecorrcoef50=mse(noisetarget50tv-noisea10tv)
testmsecorrcoef50=mse(noisetarget50tst-noisea10tst)
 
r2msetr=[trmsecorrcoef10 trmsecorrcoef20 trmsecorrcoef30 trmsecorrcoef40 trmsecorrcoef50];
r2msetv=[valmsecorrcoef10 valmsecorrcoef20 valmsecorrcoef30 valmsecorrcoef40 valmsecorrcoef50];
r2msetst=[testmsecorrcoef10 testmsecorrcoef20 testmsecorrcoef30 testmsecorrcoef40 testmsecorrcoef50];
 
figure(71)
countr2mse=[r2msetr' r2msetv' r2msetst']
ybar2mse=countr2mse(1:5,:);
bar(ybar2mse,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('MSE ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('MSE Analysis in Noisy Target Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2mse=mean(countr2mse')
%Input noise analysis based on R analysis
%%%%%%%%%% R analysis
in10tr=a10tr;
in20tr=a20tr;
in30tr=a30tr;
in40tr=a40tr;
in50tr=a50tr;
 
in10tv=a10v;
in20tv=a20v;
in30tv=a30v;
in40tv=a40v;
in50tv=a50v;
 
in10tst=a10t;
in20tst=a20t;
in30tst=a30t;
in40tst=a40t;
in50tst=a50t;
 
out10tr=tr;
out20tr=tr;
out30tr=tr;
out40tr=tr;
out50tr=tr;
 
out10tv=tav;
out20tv=tav;
out30tv=tav;
out40tv=tav;
out50tv=tav;
 
out10tst=tat;
out20tst=tat;
out30tst=tat;
out40tst=tat;
out50tst=tat;
 
trcorrcoef10=corrcoef(in10tr,out10tr)
valcorrcoef10=corrcoef(in10tv,out10tv)
testcorrcoef10=corrcoef(in10tst,out10tst)
 
noisetrcod10=diag(trcorrcoef10,1).*diag(trcorrcoef10,1)
noisevalcod10=diag(valcorrcoef10,1).*diag(valcorrcoef10,1)
noisetestcod10=diag(testcorrcoef10,1).*diag(testcorrcoef10,1)
%%%%
trcorrcoef20=corrcoef(in20tr,out20tr)
valcorrcoef20=corrcoef(in20tv,out20tv)
testcorrcoef20=corrcoef(in20tst,out20tst)
 
noisetrcod20=diag(trcorrcoef20,1).*diag(trcorrcoef20,1)
noisevalcod20=diag(valcorrcoef20,1).*diag(valcorrcoef20,1)
noisetestcod20=diag(testcorrcoef20,1).*diag(testcorrcoef20,1)
%%%%%%%%%%
trcorrcoef30=corrcoef(in30tr,out30tr)
valcorrcoef30=corrcoef(in30tv,out30tv)
testcorrcoef30=corrcoef(in30tst,out30tst)
 
noisetrcod30=diag(trcorrcoef30,1).*diag(trcorrcoef30,1)
noisevalcod30=diag(valcorrcoef30,1).*diag(valcorrcoef30,1)
noisetestcod30=diag(testcorrcoef30,1).*diag(testcorrcoef30,1)
%%%%%%%%%%%%%%%%%%%5
trcorrcoef40=corrcoef(in40tr,out40tr)
valcorrcoef40=corrcoef(in40tv,out40tv)
testcorrcoef40=corrcoef(in40tst,out40tst)
 
noisetrcod40=diag(trcorrcoef40,1).*diag(trcorrcoef40,1)
noisevalcod40=diag(valcorrcoef40,1).*diag(valcorrcoef40,1)
noisetestcod40=diag(testcorrcoef40,1).*diag(testcorrcoef40,1)
%%%%%%
trcorrcoef50=corrcoef(in50tr,out50tr)
valcorrcoef50=corrcoef(in50tv,out50tv)
testcorrcoef50=corrcoef(in50tst,out50tst)
 
noisetrcod50=diag(trcorrcoef50,1).*diag(trcorrcoef50,1)
noisevalcod50=diag(valcorrcoef50,1).*diag(valcorrcoef50,1);
noisetestcod50=diag(testcorrcoef50,1).*diag(testcorrcoef50,1)
 
r2tr=[noisetrcod10 noisetrcod20 noisetrcod30 noisetrcod40 noisetrcod50];
r2tv=[noisevalcod10 noisevalcod20 noisevalcod30 noisevalcod40 noisevalcod50];
r2tst=[noisetestcod10 noisetestcod20 noisetestcod30 noisetestcod40 noisetestcod50];
 
r2tr=sqrt(r2tr)
r2tv=sqrt(r2tv)
r2tst=sqrt(r2tst)
 
figure(72)
countr2=[r2tr' r2tv' r2tst']
ybar2=countr2(1:5,:);
bar(ybar2,'LineWidth',4,'BarWidth',1)
                
xlabel('Level of Colour Noise','FontName','Tahoma','Fontsize',24);ylabel('R ','FontName','Tahoma','Fontsize',24);
% title('Shallow Interface Depth Estimation:Test Interval','FontName','Tahoma','Fontsize',24);
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight;
title('R Analysis in Noisy Input Data');
legend('Training','Validation', 'Test','Location','northeast');
set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);axis tight
meancountr2=mean(countr2')


%----------------------------------------------------------------
 
% clc;
 clear all;
 close all;clc;
% load johdata1.dat;
% johdata=johdata1;
% %load johdata.dat;%aquifer parameter estimated from water sample and VES
% vesno=johdata(:,1);%VES no.
% BR=johdata(:,2);%Bulk resistivity(ohm-m)
% FR=johdata(:,3);%Fluid resistivity(ohm-m)
% AT=johdata(:,4);%aquifer thickness(m);
% FF=johdata(:,5);%FORMATION FACTOR;
% IFF=johdata(:,6);%Inverse of formation factor;
% %INFF=johdata(:,7);%INTRINSIC FORMATION FACTOR
% IINFF=johdata(:,7);%INVERSE OF INTRINSIC FORMATION FACTOR
% ALPHA=johdata(:,8);%alpha
% M=johdata(:,9);%m
% PORO=johdata(:,10);%POROSITY 
% HC=johdata(:,11);%HYDRAULIC CONDUCTIVITY(m/s)
% LC=johdata(:,12);%longitudinal unit conductance(S)
% TR=johdata(:,13);%TRANSVERSE RESISTANCE(TR);
% Tv=johdata(:,14);%transmissivity/t
% a=johdata(:,15);%longitude
% b=johdata(:,16);%latitude
 
%format  long;
% Generate the data set.
randn('state', 0); 
rand('state', 0); 
%ndata = 100;
%ndata=488;
ndata=211;
N=211;
% N=185;
% ndata=185;
noise = 0.1;
% a0=0.550; %initial value: /range 
% c0 = 0.06; % initial value:/ sill 
% n0=0.0015;%initial value of /nugget
% slope=0.0203;% slope value
% 
% nugget=n0;
% sill=c0;
% range=a0;
% 
% T=AT/max(AT);% Parameter of interst......................
% 
% d = variogram([a b],T,'plotit',true,'nrbins',36);
% %a0=.06;
% %c0=0.55;
% 
% [arange,ssill,nnugget,S] = variogramfit(d.distance,d.val,a0,c0,[],'solver','fminsearch','nugget',0,'plotit',true);
% 
% ssill
% arange
% nnugget
 
% load varjoh1.dat;
% x1v=varjoh1(:,1);%long
% x2v=varjoh1(:,2);%lat
% %x3v=varjoh1(:,3);%T-value
% x4v=varjoh1(:,4);%d.distance
% x5v=varjoh1(:,5);%nbins
% %x6v=varjoh1(:,6);%d.val
% %x7v=varjoh1(:,7);%range
% %x8v=varjoh1(:,8);%sill
% %x9v=varjoh1(:,9);%nugget
% x3v=T;
% x6v=d.val;
 
 
 close all;clc;
 %load WHOLE_old.dat
 load gravdeep.dat;
 %load gravshalow.dat;
 %WHOLE_old=gravshalow;
 WHOLE_old=gravdeep;
 x1v=WHOLE_old(:,1);%longitude
 x2v=WHOLE_old(:,2);%latitude
 x3v=WHOLE_old(:,3);%%BGA Anomaly
 x33v=WHOLE_old(:,5);%
%  BGA=WHOLE_old(:,3)%BGA Anomaly
%  R1A=WHOLE_old(:,4)%BGA Anomaly
%  R2A=WHOLE_old(:,5)%BGA Anomaly
%  R3A=WHOLE_old(:,6)%BGA Anomaly
%  R4A=WHOLE_old(:,7)%BGA Anomaly
%  R5A=WHOLE_old(:,8)%BGA Anomaly
%  R6A=WHOLE_old(:,9)%BGA Anomaly
%  R7A=WHOLE_old(:,10)%BGA Anomaly
%  R8A=WHOLE_old(:,11)%BGA Anomaly
%  R9A=WHOLE_old(:,12)%BGA Anomaly
%  R10A=WHOLE_old(:,13)%BGA Anomaly
%  
%n=488;%No. of Data Points
x=[x1v,x2v];
t=[x33v];
%[pn,minp,maxp]=premnmx(x);
%[tn,mint,maxt] = premnmx(t);
%xn=pn;
%tn=tn;
x=x';
t=t';
[xn,minx,maxx,tn,mint,maxt] = premnmx(x,t);
[xtrans,transMat] = prepca(xn,0.000000000000000000000000000000002);
[R,Q] = size(xtrans);
iitst = 2:4:Q;
iival = 4:4:Q;
iitr = [1:4:Q 3:4:Q];
valX = xtrans(:,iival); val.T = tn(:,iival);
testX = xtrans(:,iitst); test.T = tn(:,iitst);
xtr = xtrans(:,iitr); ttr = tn(:,iitr);
xr=xtr';
zr=ttr';
 
%net = gp(2, 'sqexp');
net = gp(2, 'ratquad');
% Initialise the parameters.
prior.pr_mean = 0;
prior.pr_var = .001;
net = gpinit(net, xr, zr, prior);
 
% Now train to find the hyperparameters.
options = foptions;
options(1) = 1;
options(14) = 300;
 
[net, options] = netopt(net, options, xr, zr, 'scg');
 
rel = exp(net.inweights);
 
 
[an, sqxr] = gpfwd(net,xr);
[y] = postmnmx(an',mint,maxt);
sqtr=2.00*sqrt(sqxr');
tr=t(:,iitr);
 
figure(26);
xis=1:106;
E1=sqtr(:,1).*ones(size(y));
%subplot(3,1,1),
errorbar(xis,y,E1,'c');grid on;axis([0 106 0 6]);hold on
plot(xis,tr,'ro-',xis, y,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 106 -6 6]);
legend('Error-bar','Target','GPreg.');
xlabel('No. of Data','Fontsize',20);
ylabel('Sediment Depth(Km)','Fontsize',20);
title('Training Interval')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',18);axis tight
% subplot(3,1,2),errorbar(xis,dm,E2,'r');grid on;axis([0 702 -1 1]);
% subplot(3,1,3),errorbar(xis,dh,E3,'b');grid on;axis([0 702 -1 1]);
% Evaluate error bars.
%[y, sigsq] = gpfwd(net, x, cninv);
%[yr, sig2,invhess] = netevfwd(mlppak(net), net, xr, zr, xr);
% [ysig] = postmnmx(sig2',mint,maxt);
% sig = sqrt(ysig);
%hess=inv(invhess);
 
 
 
%........................error analysis
dev=tr-y;
dev=dev';
dbp=dev(:,1);
msedbp=mse(dbp)
maedbp=mae(dbp)
retrain=1.0-sum((tr-y).^2)./sum(tr.^2)
dtrain=1.0-sum((tr-y).^2)./sum((y-mean(tr)+tr-mean(tr)).^2)
Rtrain=corrcoef(tr,y)
 
 
[av,sqvalx] = gpfwd(net,valX');
[avy] = postmnmx(av',mint,maxt);
sqval=2.00*sqrt(sqvalx');
[av] = postmnmx(av',mint,maxt);
tav=t(:,iival);
 
figure(27);
xis=1:52;
E2=sqtr(:,1).*ones(size(av));
 
%subplot(3,1,1),
errorbar(xis,avy,E2,'c');grid on;axis([0 52 0 6]);hold on
%title('Errobar plot');hold on
plot(xis,tav,'ro-',xis, avy,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 52 0 6]);
legend('Error-bar','Target','GPreg.');
%plot(tat);hold on;plot(ate);
xlabel('No. of Data','Fontsize',20);
ylabel('Sediment Depth(Km)','Fontsize',20);
title('Validation Interval')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',18);axis tight
 
 
devv=tav-av;
devv=devv';
vdbp=devv(:,1);
msevdbp=mse(vdbp)%mse of validation set
maevdbp=mae(vdbp)
reval=1.0-sum((tav-av).^2)./sum(tav.^2)
dval=1.0-sum((tav-av).^2)./sum((av-mean(av)+tav-mean(tav)).^2)
Rval=corrcoef(tav,av)
 
 
 
[ate,sqtestx] = gpfwd(net,testX');
[atey] = postmnmx(ate',mint,maxt);
sqtest=2.00*sqrt(sqtestx);
[ate] = postmnmx(ate',mint,maxt);
tat=t(:,iitst);
 
figure(28);
xis=1:53;
E3=sqtr(:,1).*ones(size(atey));
 
%subplot(3,1,1),
errorbar(xis,atey,E3,'c');grid on;axis([0 53 0 6]);hold on
%title('Errobar plot');hold on
plot(xis,tat,'ro-',xis, atey,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);grid on;axis([0 53 0 6]);
legend('Error-bar','Target','GPreg.');
%plot(tat);hold on;plot(ate);
xlabel('No. of Data','Fontsize',20);
ylabel('Sediment Depth(Km)','Fontsize',20);
title('Test Interval')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
devt=tat-ate;
devt=devt';
tdbp=devt(:,1);
msetdbp=mse(tdbp)%mse of test set
maetdbp=mae(tdbp)
retest=1.0-sum((tat-ate).^2)./sum(tat.^2)
dtest=1.0-sum((tat-ate).^2)./sum((tat-mean(tat)+tat-mean(tat)).^2)
Rtest=corrcoef(tat,ate)
 
 
% %disp('Press any key to continue.');
% %pause; clc;
% 
% fprintf(1, '    alpha1: %8.5f\n', net.alpha(1));
% fprintf(1, '    alpha2: %8.5f\n', net.alpha(2));
% % fprintf(1, '    alpha3: %8.5f\n', net.alpha(3));
% % fprintf(1, '    alpha4: %8.5f\n', net.alpha(4));
% % fprintf(1, '    alpha5: %8.5f\n', net.alpha(5));
% %fprintf(1, '    alpha6: %8.5f\n', net.alpha(6));
% %fprintf(1, '    alpha7: %8.5f\n', net.alpha(7));
% %fprintf(1, '    alpha8: %8.5f\n', net.alpha(8));
% %fprintf(1, '    alpha9: %8.5f\n', net.alpha(9));
% %fprintf(1, '    alpha10: %8.5f\n', net.alpha(10));
% %fprintf(1, '    alpha11: %8.5f\n', net.alpha(11));
% %fprintf(1, '    alpha12: %8.5f\n', net.alpha(12));
% 
% 
% disp('This is confirmed by looking at the corresponding weight values:')
% disp(' ');
% fprintf(1, '    %8.5f    %8.5f\n', net.w1');
% 
% lag=[1 2 ];%ph(1)-ec(2)-tds(3)-hco3(4)-cl(5)-so4(6)-no3(7)-ca(8)-mg(9)-na(10)-k(11)-rho(12)
% x1c=inv(net.alpha(1))*100;
% x2c=inv(net.alpha(2))*100;
% x3c=inv(net.alpha(3))*100;
% x4c=inv(net.alpha(4))*100;
% x5c=inv(net.alpha(5))*100;
%x6c=inv(net.alpha(6))*100;
%x7c=inv(net.alpha(7))*100;
%x8c=inv(net.alpha(8))*100;
%x9c=inv(net.alpha(9))*100;
%x10c=inv(net.alpha(10))*100;
%x11c=inv(net.alpha(11))*100;
%x12c=inv(net.alpha(12))*100;
 
% xall=[x1c x2c];
% xall=xall./max(xall)*100;
%xall=xall./20;
 
% figure(6);
% bar(lag,xall);
% set(gca,'LineWidth',3,'FontName','Malgun Gothic','Fontsize',24);axis tight;
% xlabel('Input Variables','Fontsize',24);
% ylabel('Relative contribution','Fontsize',24);
% 
% grid on;
% axis tight;
 
figure(7);
%logx1=log(net.alpha(1));
%logx2=log(net.alpha(2));
%logx3=log(net.alpha(3));
%logx4=log(net.alpha(4));
%logx5=log(net.alpha(5));
%logx6=log(net.alpha(6));
%logx7=log(net.alpha(7));
%logx8=log(net.alpha(8));
%logx9=log(net.alpha(9));
%logx10=log(net.alpha(10));
%logx11=log(net.alpha(11));
%logx12=log(net.alpha(12));
 
 
% xnew=[x1v x2v];
% 
% [n2tstp]= tramnmx(xnew',minx,maxx);
% [p2transp] = trapca(n2tstp,transMat);
% atstp = gpfwd(net,p2transp');
% [anewp] = postmnmx(atstp',mint,maxt);
% anewp=anewp';
 
 
% figure(13);
% order=1;
% plotit(x6v,anewp(:,1),order);
% MSE=mse(x6v-anewp(:,1))  % residuals
% MAE=mae(x6v-anewp(:,1))
% R=corrcoef(x6v,anewp(:,1))
% Rsquare=R.*R  
% sig1=abs(sqrt(diag(cov(anewp(:,1),S.gammahat))))
% 
% xis=1:36;
% E1=sig1(1).*ones(size(anewp(:,1)));
% %E2=sig(:,2).*ones(size(dm));
% %E3=sig(:,3).*ones(size(dh));
% 
% figure(14);
% %plot(d.distance, x6v,'o--', d.distance,anewp(:,1),'d--',d.distance,S.gammahat);
% plot(d.distance, x6v,'o--',d.distance,S.gammahat,'-rp',d.distance,anewp(:,1),'d--');
% legend( 'Experimental','Spherical','ARD-BNN');
% xlabel('Distance between observations','Fontsize',24)
% ylabel('Semivariance','Fontsize',24);
% set(gca,'LineWidth',3,'FontName','Helvetica','Fontsize',24);
% 
% MSE1=mse(x6v-S.gammahat)  % residuals
% MAE1=mae(x6v-S.gammahat)
% R1=corrcoef(x6v,S.gammahat)
% R1square=R1.*R1  
% [gamma1]=[S.gammahat];
% fit1=1.0-sum((x6v-gamma1).^2)./sum((x6v-mean(x6v)+x6v-mean(x6v)).^2)
% 
% MSE2=mse(x6v-anewp(:,1))  % residuals
% MAE2=mae(x6v-anewp(:,1))
% R2=corrcoef(x6v,anewp(:,1))
% R2square=R2.*R2  
% fit2=1.0-sum((x6v-anewp(:,1)).^2)./sum((x6v-mean(x6v)+x6v-mean(x6v)).^2)
 
%figure
figure
 
subplot(3,1,1)
%plot(tr);hold on;plot(ypred)
plot(tr,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6); hold on;plot(y,'r-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
legend('Data','GP-Pred.');
%plot(tat);hold on;plot(ate);
%xlabel('No. of Data','Fontsize',18);
ylabel('Depth[Km]','Fontsize',20);
title('Training Interval:Original Vs. GP-Predicted')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
subplot(3,1,2)
%plot(tav);hold on;plot(av);
plot(tav,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6); hold on;plot(av,'r-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
legend('Data','GP-Pred.');
%plot(tat);hold on;plot(ate);
%xlabel('No. of Data','Fontsize',18);
ylabel('Depth[Km]','Fontsize',20);
title('Validation Interval:Original Vs. GP-Predicted')
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
subplot(3,1,3)
plot(tat,'o-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6); hold on;plot(ate,'r-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6);
legend('Data','GP-Pred.');
%plot(tat);hold on;plot(ate);
xlabel('No. of Data','Fontsize',20);
ylabel('Depth[Km]','Fontsize',20);
title('Test Interval:Original Vs. GP-Predicted');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight
 
 
 
figure
plotregression(tr,y,'Regression:Training Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);
%legend('Training Interval')
xlabel('Target Depth[Km]','Fontsize',20);
ylabel('Predicted Depth[Km]','Fontsize',20);
figure
%subplot(3,2,4)
plotregression(tav,av,'Regression:Validation Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);
%legend('Validation Interval')
xlabel('Target Depth[Km]','Fontsize',20);
ylabel('Predicted Depth[Km]','Fontsize',20);
 
figure
%subplot(3,2,6)
plotregression(tat,ate,'Regression:Test Interval');
set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);
%legend('Test Interval')
xlabel('Target Depth[Km]','Fontsize',20);
ylabel('Predicted Depth[Km]','Fontsize',20);
 
load gravwhole.dat;
f=gravwhole(:,1);%longitude
g=gravwhole(:,2);%latitude
h=gravwhole(:,3);%BG value
%i=gravwhole(:,4);%altitude
figure
% Real data analysis;
vinerror10=[f,g];
[e10]= tramnmx(vinerror10',minx,maxx);
[p2e10] = trapca(e10,transMat);
[a10,SIGSQ] = gpfwd(net,p2e10');
[a10gp] = postmnmx(a10',mint,maxt);
[ysig] = postmnmx(SIGSQ',mint,maxt);
siggp =2.0*sqrt(SIGSQ);
% prediction of error
 
cninv=inv(gpcovar(net,xr));
[yp, sigsq] = gpfwd(net, p2e10', cninv);
[sigsqgp] = postmnmx(sigsq',mint,maxt);
stdgp=sqrt(sigsqgp);
 
[COVT, COVF] = gpcovar(net, p2e10');
% STDGPT=diag(sqrt(COVT));
 
xnew=[x1v x2v];
 
[n2tstp]= tramnmx(xnew',minx,maxx);
[p2transp] = trapca(n2tstp,transMat);
atstp = gpfwd(net,p2transp');
[anewp] = postmnmx(atstp',mint,maxt);
anewp=anewp';
 
 
 adata=N;
 bdata=N;
 amax=max(x1v);
 bmax=max(x2v);
 amin=min(x1v);
 bmin=min(x2v);
 ax1v = linspace(amin, amax,adata);
 bx2v = linspace(bmin, bmax,bdata);
 [ai,bi] = meshgrid(ax1v,bx2v);
 d1i = griddata(x1v,x2v,anewp,ai,bi,'v4');
% 
 % %meshz(xi,yi,zi);hold on;
% %surf(xi,yi,zi);hold on;
contourf(ai,bi,d1i);hold on;
%plot3(a,b,d,'O');
% %view(30,60);
view(0,90);
% shading interp;
% grid on;
% box on;
% axis tight;
% %axis auto;
% % axis square
 colorbar;
 set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',20);axis tight;
 xlabel('Longitude[in Degree]','Fontsize',20);
ylabel('Latitude[in degree]','Fontsize',20);
title('Spatial Variability of Sediment Depth by Using GP Model','Fontsize',20);
colormap(jet)
axis square
% hold on;
% scatter(x1v,x2v,'o');legend on;
 
%GP MODELLING FOR VARIOGRAM ANALYSIS
load WHOLE_sedi.dat;
a=WHOLE_sedi(:,1);%longitude
b=WHOLE_sedi(:,2);%latitude
c=WHOLE_sedi(:,3);%Basement Depth by svm
c=c/max(c);
d=WHOLE_sedi(:,4);%Basement Depth by gp
d=d/max(d);
e=WHOLE_sedi(:,5);%Basement Depth by ard
e=e/max(e);
f=WHOLE_sedi(:,6);%Basement Depth by anfis
f=f/max(f);
%format  long;
% Generate the data set.
randn('state', 0); 
rand('state', 0); 
%ndata = 100;
ndata=488;
noise = 0.05;
 
a0=0.550; %initial value: /range 
c0 = 0.06; % initial value:/ sill 
n0=0.0015;%initial value of /nugget
slope=0.0203;% slope value
 
nugget=n0;
sill=c0;
range=a0;
 
T=c;% Parameter of interst......................
 
dv = variogram([a b],T,'plotit',true,'nrbins',488);
 
[arange,ssill,nnugget,S] = variogramfit(dv.distance,dv.val,a0,c0,[],'solver','fminsearchbnd','nugget',0,'plotit',true);
%[a,c,n] = variogramfit(h,gammaexp,a0,c0,[],...
%                            'solver','fminsearchbnd',...
%                            'nugget',0,...
%                            'plotit',true);
ssill
arange
nnugget
 
 
x1v=a;%long
x2v=b;%lat
% %x3v=varWhole_sedi(:,3);%Basement Depth-value
% x4v=varWhole_sedi(:,4);%d.distance
% x5v=varWhole_sedi(:,5);%nbins
%x6v=varjoh1(:,6);%d.val
%x7v=varjoh1(:,7);%range
%x8v=varjoh1(:,8);%sill
%x9v=varjoh1(:,9);%nugget
x3v=T;
x4v=dv.distance;%d.distance
x5v=ones(488,1);%nbins
x6v=dv.val;%d.val
 
 
x=[x1v,x2v,x3v,x4v,x5v];
t=[x6v];
%[pn,minp,maxp]=premnmx(x);
%[tn,mint,maxt] = premnmx(t);
%xn=pn;
%tn=tn;
x=x';
t=t';
[xn,minx,maxx,tn,mint,maxt] = premnmx(x,t);
[xtrans,transMat] = prepca(xn,0.000000000000000000000000000000002);
[R,Q] = size(xtrans);
iitst = 2:4:Q;
iival = 4:4:Q;
iitr = [1:4:Q 3:4:Q];
valX = xtrans(:,iival); val.T = tn(:,iival);
testX = xtrans(:,iitst); test.T = tn(:,iitst);
xtr = xtrans(:,iitr); ttr = tn(:,iitr);
xr=xtr';
zr=ttr';
 
net = gp(5, 'sqexp');
%net = gp(5, 'ratquad');
% Initialise the parameters.
prior.pr_mean = 0;
prior.pr_var = 0.1;
net = gpinit(net, xr, zr, prior);
 
% Now train to find the hyperparameters.
options = foptions;
options(1) = 1;
options(14) = 300;
 
[net, options] = netopt(net, options, xr, zr, 'scg');
 
rel = exp(net.inweights);
 
 
an = gpfwd(net,xr);
[y] = postmnmx(an',mint,maxt);
%REAL DATA ANALYSIS
xnew=[x1v x2v x3v x4v x5v];
 
[n2tstp]= tramnmx(xnew',minx,maxx);
[p2transp] = trapca(n2tstp,transMat);
atstp = gpfwd(net,p2transp');
%atstp = svmval(p2transp',xsup,w,w0,kernel,kerneloption);
[anewp] = postmnmx(atstp',mint,maxt);
anewp=anewp';
 
gpMSE2=mse(x6v-anewp(:,1))  % residuals
gpMAE2=mae(x6v-anewp(:,1))
gpR2=corrcoef(x6v,anewp(:,1))
gpR2square=gpR2.*gpR2  
gpfit2=1.0-sum((x6v-anewp(:,1)).^2)./sum((x6v-mean(x6v)+x6v-mean(x6v)).^2)
gpfitrefit2=1.0-sum((x6v-anewp(:,1)).^2)./sum(x6v.^2)
 
% contour plot of sediment depth
a=WHOLE_sedi(:,1);%longitude
b=WHOLE_sedi(:,2);%latitude
c=WHOLE_sedi(:,3);%Basement Depth by svm
%c=c/max(c);
 
d=WHOLE_sedi(:,4);%Basement Depth by gp
%d=d/max(d);
e=WHOLE_sedi(:,5);%Basement Depth by ard
%e=e/max(e);
f=WHOLE_sedi(:,6);%Basement Depth by anfis
c=d;
%f=f/max(f);
N=488;%No. of Data Points
 
 adata=488;
 bdata=488;
 amax=max(a);
 bmax=max(b);
 amin=min(a);
 bmin=min(b);
 av = linspace(amin, amax,adata);
 bv = linspace(bmin, bmax,bdata);
 
 [X,Y] = meshgrid(av,bv);
 Z = griddata(a,b,c,X,Y,'v4');
 
 
% title('BGA[mGal]')
%
%     % calculate the sample variogram
      %v = variogram([a b],c,'plotit',false,'maxdist',488);
      v = variogram([a b],c,'plotit',true,'nrbins',488);
 
%     Fit a spherical variogram
     %v.val=anewp;
     %subplot(2,2,2)
     [dum1,dum2,dum3,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable');
     title('Variogram');axis square;
     set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',24);
     vstruct.Rs
%
%  Now use the sampled locations in a kriging
[Zhat,Zvar] = kriging(vstruct,a,b,c,X,Y);
 
figure
 
contourf(X(1,:),Y(:,1),Zhat); axis image; axis xy
shading interp;set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
colorbar;colormap(jet);
title('Sediment Depth Prediction by SVMreg')
%subplot(2,2,4)
% figure
% contourf(X,Y,Zvar); axis image
% title('MLOK Variance')
% shading interp;
% colorbar;
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
% colormap(jet);
%
%contourusing MLOK 
figure
v.val=anewp;
     %subplot(2,2,2)
     [dum1,dum2,dum3,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable');
     title('Variogram');axis square;
     set(gca,'LineWidth',3,'FontName','Tahoma','Fontsize',24);
     vstruct.Rs
%
%  Now use the sampled locations in a kriging
[Zhat1,Zvar1] = kriging(vstruct,a,b,c,X,Y);
 
figure
subplot(1,2,1)
contourf(X(1,:),Y(:,1),Zhat); %axis image; axis xy
shading interp;set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
colorbar;colormap(jet);
title('Basement Depth Predi. by SVMreg')
subplot(1,2,2)
contourf(X(1,:),Y(:,1),Zhat1); %axis image; axis xy
shading interp;set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
colorbar;colormap(jet);
title('Basement Depth Predi. by MLOK-SVMreg')
%subplot(2,2,4)
% figure
% contourf(X,Y,Zvar); axis image
% title('MLOK Variance')
% shading interp;
% colorbar;
% set(gca,'LineWidth',4,'FontName','Tahoma','Fontsize',24);
% colormap(jet);
 

