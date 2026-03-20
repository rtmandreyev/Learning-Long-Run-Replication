clear
clc

%% Load data

%load('tblResults_SlopeCoefficients_Level_SmoothL.mat')
%load('tblResults_SlopeCoefficients_Level_Break_SmoothTL.mat')
%load('tblResults_Direct_SmoothL.mat')
load('tblResults_SlopeCoefficients_Level_SmoothL.mat')
%load('tblResults_Loose_SmoothL.mat')
load('Tbl_Data.mat')
load('Tbl_Regression_Stats.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Figures

logInvGamPdf = @(x,a,b) a*log(b) - gammaln(a) + (-a-1)*log(x) - b./x;
invGamPdf = @(x,a,b) exp(logInvGamPdf(x,a,b));



%% Figure 4 (Prior Distributions)

numPlot = 10000;
xPlot1 = linspace(0,1,numPlot)';
xPlot2 = linspace(0,2,numPlot)';

fig4 = figure(4);
tiledlayout(3,1,'TileSpacing','Compact')

nexttile
plot(xPlot1,normpdf(xPlot1,rhoParams(1),sqrt(rhoParams(2))),'-k','LineWidth',3)
title('\rho')
axis tight

nexttile
plot(xPlot1,betapdf(xPlot1,gamParams(1),gamParams(2)),'-k','LineWidth',3)
title('\gamma')
axis tight

nexttile
plot(xPlot2,invGamPdf(xPlot2,sig2Params(1),sig2Params(2)),'-k','LineWidth',3)
title('\sigma^{2}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig4.WindowState = 'maximized';
pause(1)
%saveas(fig4,'tblPriors','epsc')

%% Figure G.1 (Prior Distributions for model with break)
% run this section code ONLY AFTER loading mat file: 'tblResults_SlopeCoefficients_Level_Break_SmoothTL.mat'
numPlot = 10000;
xPlot1 = linspace(0,1,numPlot)';
xPlot2 = linspace(0,2,numPlot)';

fig4 = figure(4);
tiledlayout(2,2,'TileSpacing','Compact')

nexttile
plot(xPlot1,normpdf(xPlot1,rhoParams(1),sqrt(rhoParams(2))),'-k','LineWidth',3)
title('\rho')
axis tight

nexttile
plot(xPlot1,betapdf(xPlot1,gamParams(1),gamParams(2)),'-k','LineWidth',3)
title('\gamma_1')
axis tight

nexttile
plot(xPlot1,betapdf(xPlot1,gam2Params(1),gam2Params(2)),'-k','LineWidth',3)
title('\gamma_2')
axis tight

nexttile
plot(xPlot2,invGamPdf(xPlot2,sig2Params(1),sig2Params(2)),'-k','LineWidth',3)
title('\sigma^{2}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig4.WindowState = 'maximized';
pause(1)
%saveas(fig4,'tblPriors','epsc')


%% Figure 6 (Whisker Plots) 

fig6 = figure(6);
tiledlayout(2,1,'TileSpacing','Compact')

startInd = find(~isnan(spfForecastsQ(:,1)),1);

nexttile
plot(dateVec(startInd:end),tblData(startInd:end),'-k','LineWidth',3)
hold on
for t = startInd:T
    
    plot(dateVec(t:min(T,t+4)),spfForecastsQ(t,1:min(5,5+T-t-4)),'-o',...
        'Color',[0.4 0.4 0.4],'LineWidth',2)
    
end
title('Data')

nexttile
plot(dateVec(startInd:end),tblData(startInd:end),'-k','LineWidth',3)
hold on
for t = startInd:T
    
    plot(dateVec(t:min(T,t+4)),mean(yForecastsShort(t,1:min(5,5+T-t-4),:),3,'omitnan'),...
        '-o','Color',[0.4 0.4 0.4],'LineWidth',2)
    
end
%title('Model (baseline)')
title('Model')

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig6.WindowState = 'maximized';
pause(1)
%saveas(fig6,'tblWhiskerPlots2','epsc')

%% Figure 7 (Parameter Estimates) 
fig7 = figure(7);

tiledlayout(3,1,'TileSpacing','Compact')

nexttile
plot(dateVec,mean(rhoDraws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(rhoDraws,5,2),prctile(rhoDraws,95,2)],'--k','LineWidth',3)
hold off
title('\rho')
axis tight

nexttile
plot(dateVec,mean(gamDraws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(gamDraws,5,2),prctile(gamDraws,95,2)],'--k','LineWidth',3)
hold off
title('\gamma')
axis tight

nexttile
plot(dateVec,mean(sigDraws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(sigDraws,5,2),prctile(sigDraws,95,2)],'--k','LineWidth',3)
hold off
title('\sigma')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig7.WindowState = 'maximized';
pause(1)
%saveas(fig7,'tblParameterEstimates','epsc')

%% Figure 8 (State Estimates) 
fig8 = figure(8);
%tiledlayout(2,1,'TileSpacing','Compact')

startInd = find(~isnan(muRealTimeDraws(:,1)),1);

%nexttile
plot(dateVec(startInd:end),mean(muRealTimeDraws(startInd:end,:),2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec(startInd:end),[prctile(muRealTimeDraws(startInd:end,:),5,2),...
    prctile(muRealTimeDraws(startInd:end,:),95,2)],'--k','LineWidth',3)
plot(dateVec(startInd:end),mean(muSmoothedDraws(startInd:end,:),2,'omitnan'),'Color',[0.4,0.4,0.4],'LineWidth',3)
hold off
title('\mu_{t}')
axis tight

% nexttile
% plot(dateVec(startInd:end),mean(xRealTimeDraws(startInd:end,:),2,'omitnan'),'-k','LineWidth',3)
% hold on
% plot(dateVec(startInd:end),[prctile(xRealTimeDraws(startInd:end,:),5,2),...
%     prctile(xRealTimeDraws(startInd:end,:),95,2)],'--k','LineWidth',3)
% plot(dateVec(startInd:end),mean(xSmoothedDraws(startInd:end,:),2,'omitnan'),'Color',[0.4,0.4,0.4],'LineWidth',3)
% hold off
% title('x_{t}')
% axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig8.WindowState = 'maximized';
pause(1)
%saveas(fig8,'tblStateEstimates','epsc')

%% Figure 9 (Yield Spread, 10 year - 3 month) 
% different model spec
fig14 = figure(14);

plot(dateVec,yieldData(:,40)-yieldData(:,1),'-k','LineWidth',3)
hold on
plot(dateVec(~isnan(yieldData(:,40))),yieldModel(~isnan(yieldData(:,40)),40)-...
    yieldModel(~isnan(yieldData(:,40)),1) + (mean(yieldData(:,40)-...
    yieldData(:,10)-yieldModel(:,40)+yieldModel(:,10),'omitnan')),...
    'Color',[0.4 0.4 0.4],'LineWidth',3)
hold off
legend('Data','Model')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig14.WindowState = 'maximized';
pause(1)
%saveas(fig14,'yieldSpread','epsc')

%% Figure 15 (10-year yield)

fig15 = figure(15);

plot(dateVec,yieldData(:,40),'-k','LineWidth',3)
hold on
plot(dateVec(~isnan(yieldData(:,40))),yieldModel(~isnan(yieldData(:,40)),40)+...
    (mean(yieldData(:,40)-yieldModel(:,40),'omitnan')),...
    'Color',[0.4 0.4 0.4],'LineWidth',3)
hold off
legend('Data','Model')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig15.WindowState = 'maximized';
pause(1)
%saveas(fig15,'yield10yr','epsc')

%% Figure 16 (10-year yield)

fig16 = figure(16);

fiveYearCondVar = 20*(gamDraws.*sigDraws.^2) + (1-rhoDraws.^20)./(1-rhoDraws).*((1-gamDraws).*sigDraws.^2);

plot(dateVec,mean(fiveYearCondVar,2,'omitnan'),'LineWidth',3);
hold on
plot(dateVec,[quantile(fiveYearCondVar,0.05,2),quantile(fiveYearCondVar,...
    0.95,2)],'--k','LineWidth',3);
hold off
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig16.WindowState = 'maximized';
pause(1)
%saveas(fig16,'condVar5yr','epsc')

%% Slide Figures (Posterior Slices)

ind = find(quarter(dateVec)==2 & year(dateVec) == 1981);

fig17 = figure(17);

tiledlayout(3,2)

nexttile
histogram(rhoDraws(ind,:),'Normalization','pdf')
hold on
plot(xPlot1,normpdf(xPlot1,rhoParams(1),sqrt(rhoParams(2))),'-k','LineWidth',3)
hold off
title('\rho')
axis tight

nexttile
histogram(gamDraws(ind,:),'Normalization','pdf')
hold on
plot(xPlot1,betapdf(xPlot1,gamParams(1),gamParams(2)),'-k','LineWidth',3)
hold off
title('\gamma')
axis tight

nexttile
histogram(sigDraws(ind,:).^2,'Normalization','pdf')
hold on
plot(xPlot2,invGamPdf(xPlot2,sig2Params(1),sig2Params(2)),'-k','LineWidth',3)
hold off
title('\sigma^{2}')
axis tight

nexttile
histogram(muRealTimeDraws(ind,:))
title('\mu_{t|t}')
axis tight

nexttile
histogram(xRealTimeDraws(ind,:))
title('x_{t|t}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig17.WindowState = 'maximized';
pause(1)
%saveas(fig17,'tblPosterior1981Q2','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = find(quarter(dateVec)==1 & year(dateVec) == 1995);

fig18 = figure(18);

tiledlayout(3,2)

nexttile
histogram(rhoDraws(ind,:),'Normalization','pdf')
hold on
plot(xPlot1,normpdf(xPlot1,rhoParams(1),sqrt(rhoParams(2))),'-k','LineWidth',3)
hold off
title('\rho')
axis tight

nexttile
histogram(gamDraws(ind,:),'Normalization','pdf')
hold on
plot(xPlot1,betapdf(xPlot1,gamParams(1),gamParams(2)),'-k','LineWidth',3)
hold off
title('\gamma')
axis tight

nexttile
histogram(sigDraws(ind,:).^2,'Normalization','pdf')
hold on
plot(xPlot2,invGamPdf(xPlot2,sig2Params(1),sig2Params(2)),'-k','LineWidth',3)
hold off
title('\sigma^{2}')
axis tight

nexttile
histogram(muRealTimeDraws(ind,:))
title('\mu_{t|t}')
axis tight

nexttile
histogram(xRealTimeDraws(ind,:))
title('x_{t|t}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig18.WindowState = 'maximized';
pause(1)
%saveas(fig18,'tblPosterior1995Q1','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = find(quarter(dateVec)==3 & year(dateVec) == 2007);

fig19 = figure(19);

tiledlayout(3,2)

nexttile
histogram(rhoDraws(ind,:),'Normalization','pdf')
hold on
plot(xPlot1,normpdf(xPlot1,rhoParams(1),sqrt(rhoParams(2))),'-k','LineWidth',3)
hold off
title('\rho')
axis tight

nexttile
histogram(gamDraws(ind,:),'Normalization','pdf')
hold on
plot(xPlot1,betapdf(xPlot1,gamParams(1),gamParams(2)),'-k','LineWidth',3)
hold off
title('\gamma')
axis tight

nexttile
histogram(sigDraws(ind,:).^2,'Normalization','pdf')
hold on
plot(xPlot2,invGamPdf(xPlot2,sig2Params(1),sig2Params(2)),'-k','LineWidth',3)
hold off
title('\sigma^{2}')
axis tight

nexttile
histogram(muRealTimeDraws(ind,:))
title('\mu_{t|t}')
axis tight

nexttile
histogram(xRealTimeDraws(ind,:))
title('x_{t|t}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig19.WindowState = 'maximized';
pause(1)
%saveas(fig19,'tblPosterior2007Q3','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = find(quarter(dateVec)==1 & year(dateVec) == 2019);

fig20 = figure(20);

tiledlayout(3,2)

nexttile
histogram(rhoDraws(ind,:),'Normalization','pdf')
hold on
plot(xPlot1,normpdf(xPlot1,rhoParams(1),sqrt(rhoParams(2))),'-k','LineWidth',3)
hold off
title('\rho')
axis tight

nexttile
histogram(gamDraws(ind,:),'Normalization','pdf')
hold on
plot(xPlot1,betapdf(xPlot1,gamParams(1),gamParams(2)),'-k','LineWidth',3)
hold off
title('\gamma')
axis tight

nexttile
histogram(sigDraws(ind,:).^2,'Normalization','pdf')
hold on
plot(xPlot2,invGamPdf(xPlot2,sig2Params(1),sig2Params(2)),'-k','LineWidth',3)
hold off
title('\sigma^{2}')
axis tight

nexttile
histogram(muRealTimeDraws(ind,:))
title('\mu_{t|t}')
axis tight

nexttile
histogram(xRealTimeDraws(ind,:))
title('x_{t|t}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig20.WindowState = 'maximized';
pause(1)
% saveas(fig20,'tblPosterior2019Q1','epsc')