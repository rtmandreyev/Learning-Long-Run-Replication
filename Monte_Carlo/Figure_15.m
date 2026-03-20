%% Load downward biased results without gamma updating

clear,clc;
load('Tbl_Data.mat','tblData','estimationStart','T','spfForecastsQ','yieldData')
load('Tbl_Stats_Final.mat')
data = 'New_DownwardBiasedPriors_1950_WOG';
datafile = strcat('simulationGibbs_',data,'.mat');
load(datafile);
tblSamples = muSamples + xSamples;
circFlag = 0;
estimationStart = 1; %dateVec(42) = 1961Q3, start of estimation

%% Plot Parameters

%True Params
thetaTrue = NaN(3,T);
thetaTrue(1,:) = rho*ones(T,1);
thetaTrue(2,:) = gam*ones(T,1);
thetaTrue(3,:) = sig*ones(T,1);

%Simulated Params
thetaSim = NaN(3,T,B);
thetaSim(1,:,:) = rhoMeanSim;
thetaSim(2,:,:) = gamMeanSim;
thetaSim(3,:,:) = sigMeanSim;

maxHorizon = 7;

% plot parameters, with 95% CI

close all;
figure('units','normalized','outerposition',[0 0 1 1])

dateVec = datetime(1951,1*3 + 1,1) + calquarters(0:T-1);

subplot(3,1,1)
hold on;
plot(dateVec,thetaTrue(1,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(1,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(1,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(1,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\rho')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

subplot(3,1,2)
hold on;
plot(dateVec,thetaTrue(2,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(2,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(2,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(2,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\gamma')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

subplot(3,1,3)
hold on;
plot(dateVec,thetaTrue(3,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(3,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(3,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(3,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\sigma_{\nu}')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

rhomeans_wog = nanmean(thetaSim(1,:,:),3);
%% Load normal downward biased updating case

load('Tbl_Data.mat','tblData','estimationStart','T','spfForecastsQ','yieldData')
load('Tbl_Stats_Final.mat')
data = 'New_DownwardBiasedPriors_1950';
datafile = strcat('simulationGibbs_',data,'.mat');
load(datafile);
tblSamples = muSamples + xSamples;
circFlag = 0;
estimationStart = 1; %dateVec(1) = 1951Q2

%% Plot Parameters

%True Params
thetaTrue = NaN(3,T);
thetaTrue(1,:) = rho*ones(T,1);
thetaTrue(2,:) = gam*ones(T,1);
thetaTrue(3,:) = sig*ones(T,1);

%Simulated Params
thetaSim = NaN(3,T,B);
thetaSim(1,:,:) = rhoMeanSim;
thetaSim(2,:,:) = gamMeanSim;
thetaSim(3,:,:) = sigMeanSim;

maxHorizon = 7;

% plot parameters, with 95% CI

close all;
figure('units','normalized','outerposition',[0 0 1 1])

dateVec = datetime(1951,1*3 + 1,1) + calquarters(0:T-1);

subplot(3,1,1)
hold on;
plot(dateVec,thetaTrue(1,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(1,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(1,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(1,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\rho')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

subplot(3,1,2)
hold on;
plot(dateVec,thetaTrue(2,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(2,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(2,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(2,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\gamma')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

subplot(3,1,3)
hold on;
plot(dateVec,thetaTrue(3,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(3,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(3,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(3,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\sigma_{\nu}')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

rhomeans = nanmean(thetaSim(1,:,:),3);

%% Figure 15
hold on;
plot(dateVec,thetaTrue(1,:),'Color','#808080','LineWidth',2);
plot(dateVec,rhomeans_wog,'k--','LineWidth',2)
plot(dateVec,rhomeans','k','LineWidth',2)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

legend('','Without gamma updating','With gamma updating','Location','Southeast')