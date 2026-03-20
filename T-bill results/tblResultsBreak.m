clear
clc

%% Load data

load('tblResults_SlopeCoefficients_Level_Break_SmoothTL.mat')
load('Tbl_Data.mat')
load('Tbl_Regression_Stats.mat')

rhoParams = thetOpt(1:2);
gam1Moments = thetOpt(3:4);
gam2Moments = thetOpt(5:6);
%rhoParams = [0.7422; 0.005];
%gam1Moments = [0.0235; 0.0849];
%gam2Moments = [0.0593; 0.0299];
gamFun = @(x) [(x(1)-1)/(sum(x)-2); sqrt(prod(x)/((sum(x)^2)*(sum(x)+1)))];
options = optimoptions(@fminunc,'Display','off');
if size(gam1Moments,1) < size(gam1Moments,2)
    gam1Moments = gam1Moments';
end
if size(gam2Moments,1) < size(gam2Moments,2)
    gam2Moments = gam2Moments';
end
gam1Params = fminunc(@(x) 1000*sum(((gamFun(x)-gam1Moments)./gam1Moments).^2),[2;2],options);
gam2Params = fminunc(@(x) 1000*sum(((gamFun(x)-gam2Moments)./gam2Moments).^2),[2;2],options);
sig2Params = [1.25;0.5625];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figures

logInvGamPdf = @(x,a,b) a*log(b) - gammaln(a) + (-a-1)*log(x) - b./x;
invGamPdf = @(x,a,b) exp(logInvGamPdf(x,a,b));

%% Figure 1 (SPF Whisker Plot)

fig1 = figure(1);

startInd = find(~isnan(spfForecastsQ(:,1)),1);

plot(dateVec(startInd:end),tblData(startInd:end),'-k','LineWidth',3)
hold on
for t = startInd:T
    
    plot(dateVec(t:min(T,t+4)),spfForecastsQ(t,1:min(5,5+T-t-4)),'-o',...
        'Color',[0.4 0.4 0.4],'LineWidth',2)
    
end
title('Data')

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig1.WindowState = 'maximized';
pause(1)
saveas(fig1,'tblWhiskerPlots','epsc')

%% Figure 2 (3-month T-Bill Rate)

fig2 = figure(2);

plot(dateVec,tblData,'-k','LineWidth',3)

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig2.WindowState = 'maximized';
pause(1)
saveas(fig2,'tblData','epsc')

%% Figure 4 (Prior Distributions)

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
plot(xPlot1,betapdf(xPlot1,gam1Params(1),gam1Params(2)),'-k','LineWidth',3)
title('\gamma_{1}')
axis tight

nexttile
plot(xPlot1,betapdf(xPlot1,gam2Params(1),gam2Params(2)),'-k','LineWidth',3)
title('\gamma_{2}')
axis tight

nexttile
plot(xPlot2,invGamPdf(xPlot2,sig2Params(1),sig2Params(2)),'-k','LineWidth',3)
title('\sigma^{2}')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig4.WindowState = 'maximized';
pause(1)
saveas(fig4,'tblPriorsBreak','epsc')

%% Figure 5 (Whisker Plots)

fig5 = figure(5);
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
title('Model')

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig5.WindowState = 'maximized';
pause(1)
saveas(fig5,'tblWhiskerPlots2Break','epsc')

%% Figure 6 (Parameter Estimates)

fig6 = figure(6);
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
fig6.WindowState = 'maximized';
pause(1)
saveas(fig6,'tblParameterEstimatesBreak','epsc')

%% Figure 7 (State Estimates)

fig7 = figure(7);
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
fig7.WindowState = 'maximized';
pause(1)
saveas(fig7,'tblStateEstimatesBreak','epsc')

%% Figure 14 (Yield Spread, 10-year-3 month) (FIGURE 9 MAIN PAPER)

fig14 = figure(14);

plot(dateVec,yieldData(:,40)-yieldData(:,1),'-k','LineWidth',3)
hold on
plot(dateVec(~isnan(yieldData(:,40))),(yieldModel(~isnan(yieldData(:,40)),40)-...
    yieldModel(~isnan(yieldData(:,40)),1))-...
    mean(yieldModel(:,40)-yieldModel(:,1),'omitnan')+...
    mean(yieldData(:,40)-yieldData(:,1),'omitnan'),...
    'Color',[0.4 0.4 0.4],'LineWidth',3)
hold off
legend('Data','Model')
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig14.WindowState = 'maximized';
pause(1)
saveas(fig14,'yieldSpreadBreak','epsc')

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
saveas(fig15,'yield10yrBreak','epsc')

%% Regression results for model with break
% Data cleaning

modelForecasts = squeeze(mean(yForecastsShort,3,'omitnan'));
modelForecasts(isnan(spfForecastsQ(:,1)),:) = NaN;

staggeredForecasts = spfForecastsQ;
staggeredForecastsModel = modelForecasts(:,1:5);
for ii = 2:5

    staggeredForecasts(:,ii) = circshift(staggeredForecasts(:,ii),ii-1);
    staggeredForecasts(1:ii-1,ii) = NaN;

    staggeredForecastsModel(:,ii) = circshift(staggeredForecastsModel(:,ii),ii-1);
    staggeredForecastsModel(1:ii-1,ii) = NaN;

end
forecastErrors = tblData - staggeredForecasts;
forecastRevisions = staggeredForecasts(:,1:end-1) - staggeredForecasts(:,2:end);

forecastErrorsModel = tblData - staggeredForecastsModel;
forecastRevisionsModel = staggeredForecastsModel(:,1:end-1) - staggeredForecastsModel(:,2:end);

%yieldData(isnan(yieldData(:,40)),:) = NaN;
yieldModel(isnan(yieldData(:,1:40))) = NaN;
%yieldModel(isnan(yieldData(:,40)),:) = NaN;

% endDate = 231;
% forecastErrors(endDate+1:end,:) = NaN;
% staggeredForecasts(endDate+1:end,:) = NaN;
% forecastRevisions(endDate+1:end,:) = NaN;
% yieldData(endDate+1:end,:) = NaN;
% forecastErrorsModel(endDate+1:end,:) = NaN;
% staggeredForecastsModel(endDate+1:end,:) = NaN;
% forecastRevisionsModel(endDate+1:end,:) = NaN;
% yieldModel(endDate+1:end,:) = NaN;

%% Bias

regressionResults.data.bias.Estimate = NaN(5,1);
regressionResults.data.bias.se = NaN(5,1);
regressionResults.data.bias.tStat = NaN(5,1);
regressionResults.data.bias.pValue = NaN(5,1);

regressionResults.model.bias.Estimate = NaN(5,1);
regressionResults.model.bias.se = NaN(5,1);
regressionResults.model.bias.tStat = NaN(5,1);
regressionResults.model.bias.pValue = NaN(5,1);

for ii = 1:5

    % Data

    y = forecastErrors(:,ii);
    y = y(~isnan(y));
    TCur = numel(y);
    X = ones(TCur,1);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Intercept',false,'Bandwidth',S,'Display','off');
    regressionResults.data.bias.Estimate(ii) = coeff;
    regressionResults.data.bias.se(ii) = se;
    regressionResults.data.bias.tStat(ii) = coeff./se;
    regressionResults.data.bias.pValue(ii) = 2*tcdf(-abs(coeff/se),floor(B));

    % Model

    y = forecastErrorsModel(:,ii);
    y = y(~isnan(y));
    TCur = numel(y);
    X = ones(TCur,1);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Intercept',false,'Bandwidth',S,'Display','off');
    regressionResults.model.bias.Estimate(ii) = coeff;
    regressionResults.model.bias.se(ii) = se;
    regressionResults.model.bias.tStat(ii) = coeff./se;
    regressionResults.model.bias.pValue(ii) = 2*tcdf(-abs(coeff/se),floor(B));

end

%% Autocorrelation

regressionResults.data.ar.Estimate = NaN(5,2);
regressionResults.data.ar.se = NaN(5,2);
regressionResults.data.ar.tStat = NaN(5,2);
regressionResults.data.ar.pValue = NaN(5,2);

regressionResults.model.ar.Estimate = NaN(5,2);
regressionResults.model.ar.se = NaN(5,2);
regressionResults.model.ar.tStat = NaN(5,2);
regressionResults.model.ar.pValue = NaN(5,2);

for ii = 2:5

    % Data

    y = forecastErrors(ii:end,ii);
    X = forecastErrors(1:end-ii+1,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.data.ar.Estimate(ii,:) = coeff;
    regressionResults.data.ar.se(ii,:) = se;
    regressionResults.data.ar.tStat(ii,:) = (coeff-[0;0])./se;
    regressionResults.data.ar.pValue(ii,:) = 2*tcdf(-abs((coeff-[0;0])./se),floor(B));

    % Model

    y = forecastErrorsModel(ii:end,ii);
    X = forecastErrorsModel(1:end-ii+1,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.model.ar.Estimate(ii,:) = coeff;
    regressionResults.model.ar.se(ii,:) = se;
    regressionResults.model.ar.tStat(ii,:) = (coeff-[0;0])./se;
    regressionResults.model.ar.pValue(ii,:) = 2*tcdf(-abs((coeff-[0;0])./se),floor(B));

end

%% Mincer-Zarnowitz

regressionResults.data.mz.Estimate = NaN(5,2);
regressionResults.data.mz.se = NaN(5,2);
regressionResults.data.mz.tStat = NaN(5,2);
regressionResults.data.mz.pValue = NaN(5,2);

regressionResults.model.mz.Estimate = NaN(5,2);
regressionResults.model.mz.se = NaN(5,2);
regressionResults.model.mz.tStat = NaN(5,2);
regressionResults.model.mz.pValue = NaN(5,2);

for ii = 1:5

    % Data

    y = tblData;
    X = staggeredForecasts(:,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.data.mz.Estimate(ii,:) = coeff;
    regressionResults.data.mz.se(ii,:) = se;
    regressionResults.data.mz.tStat(ii,:) = (coeff-[0;1])./se;
    regressionResults.data.mz.pValue(ii,:) = 2*tcdf(-abs((coeff-[0;1])./se),floor(B));

    % Model

    y = tblData;
    X = staggeredForecastsModel(:,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.model.mz.Estimate(ii,:) = coeff;
    regressionResults.model.mz.se(ii,:) = se;
    regressionResults.model.mz.tStat(ii,:) = (coeff-[0;1])./se;
    regressionResults.model.mz.pValue(ii,:) = 2*tcdf(-abs((coeff-[0;1])./se),floor(B));

end

%% Coibion-Gorodnichenko

regressionResults.data.cg.Estimate = NaN(5,2);
regressionResults.data.cg.se = NaN(5,2);
regressionResults.data.cg.tStat = NaN(5,2);
regressionResults.data.cg.pValue = NaN(5,2);

regressionResults.model.cg.Estimate = NaN(5,2);
regressionResults.model.cg.se = NaN(5,2);
regressionResults.model.cg.tStat = NaN(5,2);
regressionResults.model.cg.pValue = NaN(5,2);

for ii = 1:4

    % Data

    y = forecastErrors(:,ii);
    X = forecastRevisions(:,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.data.cg.Estimate(ii,:) = coeff;
    regressionResults.data.cg.se(ii,:) = se;
    regressionResults.data.cg.tStat(ii,:) = (coeff-[0;0])./se;
    regressionResults.data.cg.pValue(ii,:) = 2*tcdf(-abs((coeff-[0;0])./se),floor(B));

    % Model

    y = forecastErrorsModel(:,ii);
    X = forecastRevisionsModel(:,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.model.cg.Estimate(ii,:) = coeff;
    regressionResults.model.cg.se(ii,:) = se;
    regressionResults.model.cg.tStat(ii,:) = (coeff-[0;0])./se;
    regressionResults.model.cg.pValue(ii,:) = 2*tcdf(-abs((coeff-[0;0])./se),floor(B));

end

%% Campbell-Shiller

useInd = (1:T)';% >= 128;
nVec = [2 3 4 8 12 20 40];
numHorizons = numel(nVec);

regressionResults.data.cs.conventional.Estimate = NaN(numHorizons,1);
regressionResults.data.cs.conventional.se = NaN(numHorizons,1);
regressionResults.data.cs.conventional.tStat = NaN(numHorizons,1);
regressionResults.data.cs.conventional.pValue = NaN(numHorizons,1);

regressionResults.data.cs.contrarian.Estimate = NaN(numHorizons,1);
regressionResults.data.cs.contrarian.se = NaN(numHorizons,1);
regressionResults.data.cs.contrarian.tStat = NaN(numHorizons,1);
regressionResults.data.cs.contrarian.pValue = NaN(numHorizons,1);

regressionResults.model.cs.conventional.Estimate = NaN(numHorizons,1);
regressionResults.model.cs.conventional.se = NaN(numHorizons,1);
regressionResults.model.cs.conventional.tStat = NaN(numHorizons,1);
regressionResults.model.cs.conventional.pValue = NaN(numHorizons,1);

regressionResults.model.cs.contrarian.Estimate = NaN(numHorizons,1);
regressionResults.model.cs.contrarian.se = NaN(numHorizons,1);
regressionResults.model.cs.contrarian.tStat = NaN(numHorizons,1);
regressionResults.model.cs.contrarian.pValue = NaN(numHorizons,1);

for ii = 1:numHorizons

    % Data

    % Conventional
    
    n = nVec(ii);
    y = NaN(T,1);
    for t = 1:(T-n+1)
        y(t) = mean(tblData(t:t+n-1)) - tblData(t);
    end
    X = yieldData(:,n) - tblData;

    keepInd = (~isnan(y)&~isnan(X))&useInd;
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.data.cs.conventional.Estimate(ii) = coeff(2);
    regressionResults.data.cs.conventional.se(ii) = se(2);
    regressionResults.data.cs.conventional.tStat(ii) = (coeff(2)-1)./se(2);
    regressionResults.data.cs.conventional.pValue(ii) = ...
        2*tcdf(-abs((coeff(2)-1)./se(2)),floor(B));

    % Contrarian
    
    n = nVec(ii);
    y = NaN(T,1);
    for t = 1:(T-1)
        y(t) = yieldData(t+1,n-1) - yieldData(t,n);
    end
    X = (yieldData(:,n) - tblData)./(n-1);

    keepInd = (~isnan(y)&~isnan(X))&useInd;
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.data.cs.contrarian.Estimate(ii) = coeff(2);
    regressionResults.data.cs.contrarian.se(ii) = se(2);
    regressionResults.data.cs.contrarian.tStat(ii) = (coeff(2)-1)./se(2);
    regressionResults.data.cs.contrarian.pValue(ii) = ...
        2*tcdf(-abs((coeff(2)-1)./se(2)),floor(B));

    % Model

    % Conventional
    
    n = nVec(ii);
    y = NaN(T,1);
    for t = 1:(T-n+1)
        y(t) = mean(yieldModel(t:t+n-1,1)) - yieldModel(t,1);
    end
    X = yieldModel(:,n) - yieldModel(:,1);

    keepInd = (~isnan(y)&~isnan(X))&useInd;
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.model.cs.conventional.Estimate(ii) = coeff(2);
    regressionResults.model.cs.conventional.se(ii) = se(2);
    regressionResults.model.cs.conventional.tStat(ii) = (coeff(2)-1)./se(2);
    regressionResults.model.cs.conventional.pValue(ii) = ...
        2*tcdf(-abs((coeff(2)-1)./se(2)),floor(B));

    % Contrarian
    
    n = nVec(ii);
    y = NaN(T,1);
    for t = 1:(T-1)
        y(t) = yieldModel(t+1,n-1) - yieldModel(t,n);
    end
    X = (yieldModel(:,n) - yieldModel(:,1))./(n-1);

    keepInd = (~isnan(y)&~isnan(X))&useInd;
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    regressionResults.model.cs.contrarian.Estimate(ii) = coeff(2);
    regressionResults.model.cs.contrarian.se(ii) = se(2);
    regressionResults.model.cs.contrarian.tStat(ii) = (coeff(2)-1)./se(2);
    regressionResults.model.cs.contrarian.pValue(ii) = ...
        2*tcdf(-abs((coeff(2)-1)./se(2)),floor(B));

end





% Model

%% Model Results

disp('Model')

% Bias

[regressionResults.model.bias.Estimate(2:5)'; regressionResults.model.bias.se(2:5)';...
    regressionResults.model.bias.pValue(2:5)']

% Autocorrelation

[regressionResults.model.ar.Estimate(2:5,2)'; regressionResults.model.ar.se(2:5,2)';...
    regressionResults.model.ar.pValue(2:5,2)']

% Mincer-Zarnowitz

[regressionResults.model.mz.Estimate(2:5,2)'; regressionResults.model.mz.se(2:5,2)';...
    regressionResults.model.mz.pValue(2:5,2)']

% Coibion-Gorodnichenko

[regressionResults.model.cg.Estimate(2:5,2)'; regressionResults.model.cg.se(2:5,2)';...
    regressionResults.model.cg.pValue(2:5,2)']

% Campbell-Shiller

[regressionResults.model.cs.conventional.Estimate'; regressionResults.model.cs.conventional.se';...
    regressionResults.model.cs.conventional.pValue']

[regressionResults.model.cs.contrarian.Estimate'; regressionResults.model.cs.contrarian.se';...
    regressionResults.model.cs.contrarian.pValue']
%% Appendix G.3 and G.4, Cochrane and Piazzesi

% Data

yieldData = yieldData(:,1:20);
[T,H] = size(yieldData);
bondPrices = -(1:H).*yieldData./4;

y1 = yieldData(:,4);
bondPrices = bondPrices(:,4:4:end);
r = NaN(size(bondPrices));
f = [NaN(T,1),bondPrices(:,1:end-1) - bondPrices(:,2:end)];
for ii = 2:(H/4)
    r(:,ii) = [NaN;bondPrices(2:end,ii-1)-bondPrices(1:end-1,ii)];
end

rx = r - y1;

N = 5;
cpUnrestricted = NaN(N,6);
r2Unrestricted = NaN(N,1);
seUnrestricted = NaN(N,6);
fvUnrestricted = NaN(N,T);
fvRestricted = NaN(N,T);
for ii = 1:N

    temp = regstats2(rx(:,ii),[y1,f(:,2:5)],'linear',{'beta','hac',...
        'adjrsquare','yhat'});
    cpUnrestricted(ii,:) = temp.beta;
    r2Unrestricted(ii) = temp.adjrsquare;
    if isfield(temp.hac,'se')
        seUnrestricted(ii,:) = temp.hac.se;
    end
    fvUnrestricted(ii,:) = temp.yhat;

end

temp = regstats2(mean(rx(:,2:end),2),[y1,f(:,2:5)],'linear',{'beta',...
    'hac','yhat','adjrsquare'});
rxFitted = temp.yhat;
cpRestrictedGam = temp.beta;
r2RestrictedGam = temp.adjrsquare;
seRestrictedGam = temp.hac.se;
cpRestrictedB = NaN(N,1);
r2RestrictedB = NaN(N,1);
seRestrictedB = NaN(N,1);
for ii = 2:N
    temp = regstats2(rx(:,ii),rxFitted,'onlydata',{'beta','hac','adjrsquare','yhat'});
    cpRestrictedB(ii) = temp.beta;
    r2RestrictedB(ii) = temp.adjrsquare;
    seRestrictedB(ii) = temp.hac.se;
    fvRestricted(ii,:) = temp.yhat;
end

lineStyle = {'-o','-v','-s','-d'};
figure(2)
tiledlayout(2,1)

nexttile
hold on
for ii = 2:5
    plot(cpUnrestricted(ii,2:end),lineStyle{ii-1},'LineWidth',3)
end
hold off
legend({'2','3','4','5'})
axis tight
xticks([1 2 3 4 5])
xticklabels({'y_{1}','f_{1\rightarrow 2}','f_{2\rightarrow 3}',...
    'f_{3\rightarrow 4}','f_{4\rightarrow 5}'})
title('Unrestricted')

nexttile
hold on
for ii = 2:5
    plot(cpRestrictedB(ii)*cpRestrictedGam(2:end),lineStyle{ii-1},'LineWidth',3)
end
hold off
legend({'2','3','4','5'})
axis tight
xticks([1 2 3 4 5])
xticklabels({'y_{1}','f_{1\rightarrow 2}','f_{2\rightarrow 3}',...
    'f_{3\rightarrow 4}','f_{4\rightarrow 5}'})
title('Restricted')

set(findall(gcf,'-property','FontSize'),'FontSize',20)

% Model

yieldModel = yieldModel(:,1:20);
[T,H] = size(yieldModel);
meanInd = (max(sum(isnan(yieldData)))+1:T);
yieldModel = (yieldModel-mean(yieldModel,'omitnan')) + ...
    mean(yieldData(meanInd,1:H),'omitnan');
bondPrices = -(1:H).*yieldModel./4;

y1 = yieldModel(:,4);
bondPrices = bondPrices(:,4:4:end);
r = NaN(size(bondPrices));
f = [NaN(T,1),bondPrices(:,1:end-1) - bondPrices(:,2:end)];
for ii = 2:(H/4)
    r(:,ii) = [NaN;bondPrices(2:end,ii-1)-bondPrices(1:end-1,ii)];
end

rx = r - y1;

N = 5;
cpUnrestrictedModel = NaN(N,6);
r2UnrestrictedModel = NaN(N,1);
seUnrestrictedModel = NaN(N,6);
fvUnrestrictedModel = NaN(N,T);
fvUnrestrictedModelD = NaN(N,T);
fvRestrictedModel = NaN(N,T);
fvRestrictedModelD = NaN(N,T);
for ii = 1:N

    temp = regstats2(rx(:,ii),[y1,f(:,2:5)],'linear',{'beta','hac',...
        'adjrsquare','yhat'});
    cpUnrestrictedModel(ii,:) = temp.beta;
    r2UnrestrictedModel(ii) = temp.adjrsquare;
    if isfield(temp.hac,'se')
        seUnrestrictedModel(ii,:) = temp.hac.se;
    end
    fvUnrestrictedModel(ii,:) = temp.yhat;
    fvUnrestrictedModelD(ii,:) = cpUnrestricted(ii,:)*[ones(size(y1)),y1,f(:,2:5)]';

end

temp = regstats2(mean(rx(:,2:end),2),[y1,f(:,2:5)],'linear',{'beta',...
    'hac','yhat','adjrsquare'});
rxFitted = temp.yhat;
cpRestrictedGamModel = temp.beta;
r2RestrictedGamModel = temp.adjrsquare;
seRestrictedGamModel = temp.hac.se;
cpRestrictedBModel = NaN(N,1);
r2RestrictedBModel = NaN(N,1);
seRestrictedBModel = NaN(N,1);
for ii = 2:N
    temp = regstats2(rx(:,ii),rxFitted,'onlydata',{'beta','hac','adjrsquare','yhat'});
    cpRestrictedBModel(ii) = temp.beta;
    r2RestrictedBModel(ii) = temp.adjrsquare;
    seRestrictedBModel(ii) = temp.hac.se;
    fvRestrictedModel(ii,:) = temp.yhat;
    fvRestrictedModelD(ii,:) = cpRestrictedB(ii)*(cpRestrictedGam'*...
        [ones(size(y1)),y1,f(:,2:5)]');
end

figure(3)
tiledlayout(2,1)

nexttile
hold on
for ii = 2:5
    plot(cpUnrestrictedModel(ii,2:end),lineStyle{ii-1},'LineWidth',3)
end
hold off
legend({'2','3','4','5'})
axis tight
xticks([1 2 3 4 5])
xticklabels({'y_{1}','f_{1\rightarrow 2}','f_{2\rightarrow 3}',...
    'f_{3\rightarrow 4}','f_{4\rightarrow 5}'})
title('Unrestricted')

nexttile
hold on
for ii = 2:5
    plot(cpRestrictedBModel(ii)*cpRestrictedGamModel(2:end),lineStyle{ii-1},'LineWidth',3)
end
hold off
legend({'2','3','4','5'})
axis tight
xticks([1 2 3 4 5])
xticklabels({'y_{1}','f_{1\rightarrow 2}','f_{2\rightarrow 3}',...
    'f_{3\rightarrow 4}','f_{4\rightarrow 5}'})
title('Restricted')

set(findall(gcf,'-property','FontSize'),'FontSize',20)

save('regressionResults')

