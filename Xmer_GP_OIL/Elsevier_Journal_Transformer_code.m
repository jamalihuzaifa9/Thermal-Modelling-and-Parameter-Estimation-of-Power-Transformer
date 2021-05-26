clc, close all,clear all;
temp_amb=xlsread('MSEDCLdata.xlsx','J2:J46612'); % given ambient temperature data
OTS = xlsread('MSEDCLdata.xlsx','L2:L46612');  %oil temp of streaming data
load = xlsread('MSEDCLdata.xlsx','Q2:Q46612');  % given loading
WTS = xlsread('MSEDCLdata.xlsx','K2:K46612');

N = 5000;
xdata=[load(2:N,1),temp_amb(2:N,1),OTS(1:N-1,1)];
y = OTS(2:N,1);
gpMdl = fitrgp(xdata,y,'FitMethod','sr','PredictMethod','exact','BasisFunction','linear','KernelFunction',...
    'ardsquaredexponential','ComputationMethod','qr','Standardize',1);
Kernel_properties = gpMdl.KernelInformation;
tini = [OTS(1,1)];
yint = [];
M=20000;
for i = 1:M
      %m_pu_gp(i,1) = ((0.13573*10^(-5))*exp((2797.3)/(tini(i,1)+273)))/((0.13573*10^(-5))*exp((2797.3)/(48+273)));
     [est(i,1),~,yi] =  predict(gpMdl,[load(i),temp_amb(i),tini(i,1)]);
     tini(i+1,1) = est(end,1);%[tini;est(end,1)];
     yint=[yint;yi];
 %xdata = [xdata;[load_test(i),temp_amb(i),tini(i,1)]];
 %y = [y;est(end,1)];
 %gpMdl = fitrgp(xdata,y);
 end
 
 t = linspace(1,M,M)';
 figure(2)
 hold on
 plot(OTS(1:M),'--','LineWidth',2)
 plot(est,'LineWidth',2)
 patch([t;flipud(t)],[yint(:,1);flipud(yint(:,2))],'k','FaceAlpha',0.1);
 legend({'oil temp test','GP_{est}','confidence interval'})