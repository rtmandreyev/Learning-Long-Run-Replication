
clear
clc
load('Tbl_Data.mat')


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


%% Figure 2 (3-month T-Bill Rate)

fig2 = figure(2);

plot(dateVec,tblData,'-k','LineWidth',3)

set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig2.WindowState = 'maximized';
pause(1)

%%


load('Gdp_Data.mat')

%% Figure 3
dateVecAnnual = (1975:2019)';
TAnnual = numel(dateVecAnnual);
%meanForecasts  = squeeze(mean(yGrowthForecastsADraws,3,'omitnan'));
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


fig3 = figure(3);


startInd = find(~isnan(cboData(:,1)),1);
plotRangeY = [min(min(min(cboPlot))),max(max(max(cboPlot)))];


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



set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig3.WindowState = 'maximized';
pause(1)
