clear
clc

load('Tbl_Data.mat')
load('Tbl_Regression_Stats.mat')
%load('tblResults_SlopeCoefficients_Level_SmoothL.mat','yForecastsShort',...
%     'yieldModel')
%load('tblResults_SlopeCoefficients_Level_Break_SmoothTL.mat','yForecastsShort',...
%     'yieldModel')
%load('tblResults_Direct_SmoothL.mat','yForecastsShort',...
%    'yieldModel')
load('tblResults_SlopeCoefficients_Level_SmoothL.mat','yForecastsShort',...
   'yieldModel')
%load('tblResults_Loose_SmoothL.mat','yForecastsShort',...
%    'yieldModel')

%% Data cleaning

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




%% Data Results

disp('Data')

% Bias

[regressionResults.data.bias.Estimate(2:5)'; regressionResults.data.bias.se(2:5)';...
    regressionResults.data.bias.pValue(2:5)']

% Autocorrelation

[regressionResults.data.ar.Estimate(2:5,2)'; regressionResults.data.ar.se(2:5,2)';...
    regressionResults.data.ar.pValue(2:5,2)']

% Mincer-Zarnowitz

[regressionResults.data.mz.Estimate(2:5,2)'; regressionResults.data.mz.se(2:5,2)';...
    regressionResults.data.mz.pValue(2:5,2)']

% Coibion-Gorodnichenko

[regressionResults.data.cg.Estimate(2:5,2)'; regressionResults.data.cg.se(2:5,2)';...
    regressionResults.data.cg.pValue(2:5,2)']

% Campbell-Shiller

[regressionResults.data.cs.conventional.Estimate'; regressionResults.data.cs.conventional.se';...
    regressionResults.data.cs.conventional.pValue']

[regressionResults.data.cs.contrarian.Estimate'; regressionResults.data.cs.contrarian.se';...
    regressionResults.data.cs.contrarian.pValue']


% add comments
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


%% Cochrane and Piazzesi

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