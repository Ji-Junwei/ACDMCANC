clc;clear;close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');

%% configuration
Fs = 16000; % sampling frequency
T  = 120;     % time
t  = 0:1/Fs:T;
N  = length(t);

load("simulation path/SecondaryPath_6x6.mat");
load("simulation path/PrimaryPath_1x6.mat");

PrimaryPath = Primary_path;
SecondaryPath = Secondary_path;


%% system parameters
wLen = 512;  % local control filter length
sLen = 256;  % secondary path length
Numnode = 6; % number of node
cLen = 33;  % compensate filter length
muw = 1e-6; % control filter step size
muc = 1e-5; % compensate filter step size

%% noise generation
noise = randn(N,1);  % random noise
low = 100;
high = 1000;
fil = fir1(63,[2*low/Fs 2*high/Fs]);
Ref = filter(fil,1,noise);   % reference

for i = 1:Numnode
  Dis(i,:) = filter(PrimaryPath(i,:),1,Ref);   % Disturbance         
end

Ref = awgn(Ref,40,'measured');

%% centralized control
CMcANC = McANC_FxLMS_SIMO(wLen,SecondaryPath,sLen,Numnode,Numnode,N,Dis);
[e_CMANC,CMcANC] = McFxLMS_SIMO_166(CMcANC,Ref,muw);

%% MGDFxLMS
MGDFxLMS = DMANC_CompensateSP(wLen,SecondaryPath,sLen,Numnode,N,Dis,cLen);
[~,MGDFxLMS] = CompensateSP(MGDFxLMS,muc);
[e_MGDFxLMS,MGDFxLMS] = DMANC_gradient_166(MGDFxLMS,Ref,muw);

%% proposed method
Wcsubopt = zeros(Numnode,(wLen+cLen-1));
DMCANC = AC_DMCANC(wLen,SecondaryPath,sLen,Numnode,N,Dis,Ref,cLen,Wcsubopt);
[~,DMCANC] = CompensateSecP(DMCANC,muc);

alpha = [800 800 800 800 800 800];
muw = 1e-6;
DMCANC_AC = DMCANC;
DMCANC_SC = DMCANC;
[e_SCDMCANC,Iscomm_SC,DMCANC_SC] = SC_DMCANC_166(DMCANC_SC,muw,alpha,Wcsubopt);
[e_ACDMCANC,Iscomm,DMCANC_AC] = AC_DMCANC_166(DMCANC_AC,muw,alpha,Wcsubopt);

%% plot figure

figure;
for i = 1:6
    dis22 = smooth((Dis(i,1:T*Fs).^2),2000);
    ecmanc_adaptive = smooth((e_CMANC(i,1:T*Fs).^2),2000);
    emgdfxlms = smooth((e_MGDFxLMS(i,1:T*Fs).^2),2000);
    eSCDMCANC = smooth((e_SCDMCANC(i,1:T*Fs).^2),2000);
    eACDMCANC = smooth((e_ACDMCANC(i,1:T*Fs).^2),2000);

    mse_adaptive = 10*log10(ecmanc_adaptive./dis22);
    mse_mgdfxlms = 10*log10(emgdfxlms./dis22);
    mse_SCDMCANC = 10*log10(eSCDMCANC./dis22);
    mse_ACDMCANC = 10*log10(eACDMCANC./dis22);

    subplot(3,2,i);
    plot(smooth(mse_adaptive(100:end-1000,1),5000));
    hold on;
    plot(smooth(mse_mgdfxlms(100:end-1000,1),5000));
    hold on;
    plot(smooth(mse_SCDMCANC(100:end-1000,1),5000));
    hold on;
    plot(smooth(mse_ACDMCANC(100:end-1000,1),5000));
    legend('Centralized','MGDFxLMS','SCDMCANC','ACDMCANC');
    % axis([0 inf -inf 10]);
    grid on;
end

nse_adaptive = zeros(Numnode,T*Fs);
nse_mgdfxlms = zeros(Numnode,T*Fs);
nse_SCDMCANC = zeros(Numnode,T*Fs);
nse_ACDMCANC = zeros(Numnode,T*Fs);

for i = 1:6
    dis22 = smooth((Dis(i,1:T*Fs).^2),2000);
    ecmanc_adaptive = smooth((e_CMANC(i,1:T*Fs).^2),2000);
    emgdfxlms = smooth((e_MGDFxLMS(i,1:T*Fs).^2),2000);
    eSCDMCANC = smooth((e_SCDMCANC(i,1:T*Fs).^2),2000);
    eACDMCANC = smooth((e_ACDMCANC(i,1:T*Fs).^2),2000);

    nse_adaptive(i,:)    = 10*log10(ecmanc_adaptive./dis22);
    nse_mgdfxlms(i,:)    = 10*log10(emgdfxlms./dis22); 
    nse_SCDMCANC(i,:) = 10*log10(eSCDMCANC./dis22);   
    nse_ACDMCANC(i,:) = 10*log10(eACDMCANC./dis22);    
end

mse_adaptive = mean(nse_adaptive,1);
mse_mgdfxlms = mean(nse_mgdfxlms,1);
mse_SCDMCANC = mean(nse_SCDMCANC,1);
mse_ACDMCANC = mean(nse_ACDMCANC,1);

figure;
plot(smooth(mse_adaptive(100:end-1000),5000));
hold on;
plot(smooth(mse_mgdfxlms(100:end-1000),5000));
hold on;
plot(smooth(mse_SCDMCANC(100:end-1000),5000));
hold on;
plot(smooth(mse_ACDMCANC(100:end-1000),5000));
legend('Centralized','MGDFxLMS','SCDMCANC','ACDMCANC');
% axis([0 inf -inf 10]);
grid on;

figure;
for j = 1:6
    subplot(3,2,j);
    plot(Iscomm(j,:));
end

%save('case','MGDFxLMS',"CMcANC","DMCANC_SC",'e_MGDFxLMS','e_CMANC','e_SCDMCANC','Ref','Dis',"alpha","DMCANC_AC",'e_ACDMCANC','Iscomm','Iscomm_SC');