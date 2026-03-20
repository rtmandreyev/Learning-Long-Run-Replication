clear
clc

%% Load data

load('gdpResults_AllCoefficients_Levels5.mat')
load('Gdp_Data.mat')
load('Gdp_Regression_Stats.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figures

logInvGamPdf = @(x,a,b) a*log(b) - gammaln(a) + (-a-1)*log(x) - b./x;
invGamPdf = @(x,a,b) exp(logInvGamPdf(x,a,b));

%% Figure 10 (Prior Distributions) 

numPlot = 10000;
xPlot1 = linspace(0,1,numPlot)';
xPlot2 = linspace(-1,0,numPlot)';
xPlot3 = linspace(-0.03,0.05,numPlot)';
xPlot4 = linspace(0,2e-3,numPlot)';

fig11 = figure(10);
tiledlayout(3,2,'TileSpacing','Compact')

nexttile
plot(xPlot1,normpdf(xPlot1,rhoParams(1),sqrt(rhoParams(2))),'-k','LineWidth',3)
title('\rho_{1}+\rho_{2}')
axis tight

nexttile
plot(xPlot2,normpdf(xPlot2,rho2Params(1),sqrt(rho2Params(2))),'-k','LineWidth',3)
title('\rho_{2}')
axis tight

nexttile
plot(xPlot1,betapdf(xPlot1,gamParams(1),gamParams(2)),'-k','LineWidth',3)
title('\gamma')
axis tight

nexttile
plot(xPlot3,normpdf(xPlot3,muParams(1),sqrt(muParams(2))),'-k','LineWidth',3)
title('\mu')
axis tight

nexttile
plot(100*xPlot4,invGamPdf(xPlot4,sig2Params(1),sig2Params(2)),'-k','LineWidth',3)
title('100\times\sigma^{2}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig11.WindowState = 'maximized';
pause(1)
saveas(fig11,'gdpPriors5','epsc')

%% Figure 12 (Parameter Estimates) 
fig12 = figure(12);
tiledlayout(2,3,'TileSpacing','Compact')

nexttile
plot(dateVec,mean(rho1Draws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(rho1Draws,5,2),prctile(rho1Draws,95,2)],'--k','LineWidth',3)
hold off
title('\rho_{1}')
axis tight

nexttile
plot(dateVec,mean(rho2Draws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(rho2Draws,5,2),prctile(rho2Draws,95,2)],'--k','LineWidth',3)
hold off
title('\rho_{2}')
axis tight

nexttile
plot(dateVec,mean(gamDraws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(gamDraws,5,2),prctile(gamDraws,95,2)],'--k','LineWidth',3)
hold off
title('\gamma')
axis tight

nexttile
plot(dateVec,mean(muDraws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(muDraws,5,2),prctile(muDraws,95,2)],'--k','LineWidth',3)
hold off
title('\mu')
axis tight

nexttile
plot(dateVec,mean(sigDraws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(sigDraws,5,2),prctile(sigDraws,95,2)],'--k','LineWidth',3)
hold off
title('\sigma')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig12.WindowState = 'maximized';
pause(1)
saveas(fig12,'gdpParameterEstimates5','epsc')

%% Figure 11 (Whisker Plots) 
dateVecAnnual = (1975:2019)';
TAnnual = numel(dateVecAnnual);
meanForecasts  = squeeze(mean(yGrowthForecastsADraws,3,'omitnan'));
gdpInitial = NaN(size(dateVecAnnual));
T = size(gdpData,1);
curInd = 1;
curVintage = 42;
for t = 116:4:T
    gdpInitial(curInd) = 100*(mean(gdpData(t-3:t,curVintage))./...
        mean(gdpData(t-7:t-4,curVintage))-1);
    curInd = curInd+1;
    curVintage = curVintage + 4;
end

cboPlot = [gdpInitial,cboData];
modelPlot = [gdpInitial,meanForecasts];

fig11 = figure(11);
tiledlayout(2,1,'TileSpacing','Compact')

startInd = find(~isnan(cboData(:,1)),1);
plotRangeY = [min(min(min(cboPlot)),min(min(modelPlot))),...
    max(max(max(cboPlot)),max(max(modelPlot)))];

nexttile
plot(dateVecAnnual,[NaN;forecastEvaluationData(1:end-1)],'-k','LineWidth',3)
hold on
plot(dateVecAnnual,cboPlot(:,1),'--k','LineWidth',3)
legend('Most recent vintage','Initial release','Location','Best','AutoUpdate','off')
for t = 1:(TAnnual-1)
    
    plot(dateVecAnnual(t:min(TAnnual,t+6)),cboPlot(t,1:min(7,7+TAnnual-t-6)),'-o',...
        'Color',[0.4 0.4 0.4],'LineWidth',2)
    
end
title('Data')
axis tight
ylim([plotRangeY(1),plotRangeY(2)])

nexttile
plot(dateVecAnnual,[NaN;forecastEvaluationData(1:end-1)],'-k','LineWidth',3)
hold on
plot(dateVecAnnual,cboPlot(:,1),'--k','LineWidth',3)
for t = 1:(TAnnual-1)
    
    plot(dateVecAnnual(t:min(TAnnual,t+6)),modelPlot(t,1:min(7,7+TAnnual-t-6)),'-o',...
        'Color',[0.4 0.4 0.4],'LineWidth',2)
    
end
title('Model')
axis tight
ylim([plotRangeY(1),plotRangeY(2)])

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig11.WindowState = 'maximized';
pause(1)
saveas(fig11,'gdpWhiskerPlots5','epsc')

%% Figure (Cycle Estimate) (not there in the paper)

fig11 = figure(11);
%tiledlayout(2,3,'TileSpacing','Compact')

plot(dateVec,mean(xRealTimeDraws,2,'omitnan'),'-k','LineWidth',3)
hold on
plot(dateVec,[prctile(xRealTimeDraws,5,2),prctile(xRealTimeDraws,95,2)],'--k','LineWidth',3)
plot(dateVec,[NaN(estimationStart-1,1);mean(xSmoothedDraws(estimationStart:end,:),2,'omitnan')],...
   'Color',[0.4 0.4 0.4],'LineWidth',3)
% plot(dateVec(estimationStart:end),mean(xRealTimeDraws(estimationStart:end,:),2,'omitnan'),'-k','LineWidth',3)
% hold on
% plot(dateVec(estimationStart:end),[prctile(xRealTimeDraws(estimationStart:end,:),5,2),prctile(xRealTimeDraws(estimationStart:end,:),95,2)],'--k','LineWidth',3)
% plot(dateVec(estimationStart:end),mean(xSmoothedDraws(116-51+1:end,:),2,'omitnan'),...
%     'Color',[0.4 0.4 0.4],'LineWidth',3)
hold off
title('x_{t}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig11.WindowState = 'maximized';
pause(1)
saveas(fig11,'gdpCycleEstimates5','epsc')