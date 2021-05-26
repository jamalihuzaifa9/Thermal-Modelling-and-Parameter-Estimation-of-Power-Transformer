
clear all
close all
clc 
rng('default')

%delta_theta_oil_plot=xlsread('PARAMETER_16_MARCH_2021.xlsx','B1:B26'); 
%Toil_plot = xlsread('PARAMETER_16_MARCH_2021.xlsx','A1:A26');  
%n_plot = xlsread('PARAMETER_16_MARCH_2021.xlsx','C1:C25'); 
%delta_theta_oil_GP=xlsread('PARAMETER_16_MARCH_2021.xlsx','E1:E26'); 
%Toil_GP = xlsread('PARAMETER_16_MARCH_2021.xlsx','D1:D26');  
%n_GP = xlsread('PARAMETER_16_MARCH_2021.xlsx','F1:F26');  

N=6000;
temp_amb=27+0.01*randn(N,1); %
OTS = 26.5; 
load1=0.4+0.01*rand(N/3,1);
load2=0.4+0.01*rand(N/3,1);
load3=0.4+0.01*rand(N/3,1);
load = [load1;load2;load3];

M=6000;
load4=0.7+0.01*rand(M/3,1);
load5=0.5+0.01*rand(M/3,1);
load6=0.2+0.01*rand(M/3,1);
load_test = [load5;load4;load6];
temp_amb_test = 26+0.01*randn(M,1);

load7=0.4+0.05*rand(M/3,1);
load8=0.4+0.05*rand(M/3,1);
load9=0.4+0.05*rand(M/3,1);
load_test_1 = [load7;load8;load9];
temp_amb_test1 = 30+0.01*randn(M,1);

Mu = ((0.13573*10^(-5))*exp((2797.3)./(OTS+273)))./((0.13573*10^(-5))*exp((2797.3)/(70+273)));
R = 37000/5400; 
A = (1+(R*(load.^2)))/(1+R); 
A2 = (1+(R*(load_test.^2)))/(1+R); 
A3 = (1+(R*(load_test_1.^2)))/(1+R); 
%% Calculating OIL temperature
t_initial = 53.65;
n = 0.25;
Toil= 100; 
delta_theta_oil = 70;
for i=1:N
    m_pu(i,1) = ((0.13573*10^(-5))*exp((2797.3)/(t_initial+273)))/((0.13573*10^(-5))*exp((2797.3)/(70+273)));
    oil_temp_calculated(i,1)=((5*(A(i,1).*delta_theta_oil ./Toil))-(((t_initial-temp_amb(i,1)).^(1+n)/(Toil*delta_theta_oil .^n.*m_pu(i,1).^n))*5))+t_initial; 
    oil_temp_calculated = abs(oil_temp_calculated);
    t_initial=oil_temp_calculated(end);
end
oil_temp_calculated=awgn(oil_temp_calculated,10e16,0);
figure(1)
plot(oil_temp_calculated,'linewidth',2)
legend('oil_{train}')

t_initial_test = 26; 
n = 0.25;
Toil= 100;
delta_theta_oil = 70;
for i=1:M
    m_pu_test(i,1) = ((0.13573*10^(-5))*exp((2797.3)/(t_initial_test+273)))/((0.13573*10^(-5))*exp((2797.3)/(70+273)));
    oil_temp_test(i,1)=((5*(A2(i,1).*delta_theta_oil ./Toil))-(((t_initial_test-temp_amb_test(i,1)).^(1+n)/(Toil*delta_theta_oil .^n.*m_pu_test(i,1).^n))*5))+t_initial_test; % 5 is the sampling time in mins
    oil_temp_test = abs(oil_temp_test);
    t_initial_test=oil_temp_test(end);
end


T = @(x)((5*(A(1:N-1,1).*x(2)./x(1)))-((((oil_temp_calculated(1:N-1))-temp_amb(2:N)).^(1+x(3))./...
    (x(2).^(x(3)).*x(1).*m_pu(1:N-1,1).^x(3)))*5))+(oil_temp_calculated(1:N-1))-...
    (oil_temp_calculated(2:N,1));
x0 = [10,700,1];
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt'); 
x = lsqnonlin(T,x0,[],[],options); 
n1 = abs(x(3));
Toil1=abs(x(1));
delta_theta_oil1 = abs(x(2)); 

xdata=[load(2:N,1),temp_amb(2:N,1),oil_temp_calculated(1:N-1,1)];
y = oil_temp_calculated(2:N,1);
gpMdl = fitrgp(xdata,y,'FitMethod','sr','PredictMethod','sr','BasisFunction','pureQuadratic','KernelFunction',...
    'squaredexponential','ComputationMethod','qr','Standardize',1);
  
tini = 26;
 
yint1=[];
 
for i = 1:M
      m_pu_gp(i,1) = ((0.13573*10^(-5))*exp((2797.3)/(tini(i,1)+273)))/((0.13573*10^(-5))*exp((2797.3)/(48+273)));
     [est(i,1),~,yi] =  predict(gpMdl,[load_test(i),temp_amb_test(i),tini(i,1)]);
     tini(i+1,1) = est(end,1);%[tini;est(end,1)];
     yint1=[yint1;yi];
 %xdata = [xdata;[load_test(i),temp_amb(i),tini(i,1)]];
 %y = [y;est(end,1)];
 %gpMdl = fitrgp(xdata,y);
 end
 
 t = linspace(1,M,M)';
 figure(2)
 hold on
 plot(oil_temp_test,'--','LineWidth',2)
 plot(est,'LineWidth',2)
 patch([t;flipud(t)],[yint1(:,1);flipud(yint1(:,2))],'k','FaceAlpha',0.1);
 legend({'oil_temp_test','GP_{est}','confidence interval'})
 
 T1 = @(x)((5*(A2(1:M-1,1).*x(2)./x(1)))-((((est(1:M-1))-temp_amb_test(2:M)).^(1+x(3))./...
     (x(2).^(x(3)).*x(1).*m_pu_gp(1:M-1,1).^x(3)))*5))+(est(1:M-1))-...
     (est(2:M,1));
 x0 = [10,700,1];  
 options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt'); %selection of Optimisation method
 x1 = lsqnonlin(T1,x0,[],[],options); % parameter estimation using curve fitting
 n2 = abs(x1(3));
 Toil2=abs(x1(1)); 
 delta_theta_oil2 = abs(x1(2));
 
% figure(3)
% hold on
% plot(delta_theta_oil_plot)
% plot(delta_theta_oil_GP)
% plot(Toil_plot)
% plot(Toil_GP)
% plot(Toil*ones(26),'--')
% plot(delta_theta_oil*ones(26),'--')
% legend({'\Delta oil_{train}','\Delta oil_{GP}','\tau_{train}','\tau_{GP}','\tau_{ref}','\Delta oil_{ref}'})
% 
% figure(4)
% hold on
% plot(n_plot)
% plot(n_GP)
% plot(n*ones(26),'--')
% legend({'n_{train}','n_{GP}','n_{ref}'})