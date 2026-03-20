clear
clc

%% T-bill model

load('Tbl_Data.mat')

staggeredForecasts = spfForecastsQ;

for ii = 2:5

    staggeredForecasts(:,ii) = circshift(staggeredForecasts(:,ii),ii-1);
    staggeredForecasts(1:ii-1,ii) = NaN;

    

end
forecastErrors = tblData - staggeredForecasts;
forecastRevisions = staggeredForecasts(:,1:end-1) - staggeredForecasts(:,2:end);
%%

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

    

end

%% Autocorrelation

regressionResults.data.ar.Estimate = NaN(5,2);
regressionResults.data.ar.se = NaN(5,2);
regressionResults.data.ar.tStat = NaN(5,2);
regressionResults.data.ar.pValue = NaN(5,2);



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

   

end

%% Mincer-Zarnowitz

regressionResults.data.mz.Estimate = NaN(5,2);
regressionResults.data.mz.se = NaN(5,2);
regressionResults.data.mz.tStat = NaN(5,2);
regressionResults.data.mz.pValue = NaN(5,2);



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

    
    

end

%% Coibion-Gorodnichenko

regressionResults.data.cg.Estimate = NaN(5,2);
regressionResults.data.cg.se = NaN(5,2);
regressionResults.data.cg.tStat = NaN(5,2);
regressionResults.data.cg.pValue = NaN(5,2);


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

    
end

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

%% GDP regressions

load('Gdp_Data.mat')

staggeredCBO = cboData;
maxHorizon = 10;
for ii = 0:maxHorizon
 
    staggeredCBO(:,ii+1) = circshift(cboData(:,ii+1),ii);
    staggeredCBO(1:ii,ii+1) = NaN;

end

forecastEvaluationData(end) = NaN;
cboForecastErrors = forecastEvaluationData - staggeredCBO;

T = size(cboData,1);

%% Compute forecast anomalies

biasMat = NaN(maxHorizon+1,1);
biasMat_se = biasMat;
biasMat_p = biasMat;

arMat = NaN(maxHorizon+1,2);
arMat_se = arMat;
arMat_p = arMat;

mzMat = NaN(maxHorizon+1,2);
mzMat_se = mzMat;
mzMat_p = mzMat;

cgMat = NaN(maxHorizon+1,2);
cgMat_se = cgMat;
cgMat_p = cgMat;

for ii = 1:(maxHorizon+1)

    % Bias

    y = cboForecastErrors(ii:end,ii);
    y = y(~isnan(y));
    TCur = size(y,1);
    X = ones(TCur,1);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Intercept',false,'Bandwidth',S,'Display','off');
    biasMat(ii,1) = coeff;
    biasMat_se(ii,1) = se;
    biasMat_p(ii,1) = 2*tcdf(-abs(coeff./se),floor(B));

    % Autocorrelation of forecast errors

    y = cboForecastErrors(ii+1:end,ii);
    X = cboForecastErrors(1:end-ii,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    arMat(ii,:) = coeff;
    arMat_se(ii,:) = se;
    arMat_p(ii,:) = 2*tcdf(-abs(coeff./se),floor(B));

    % Mincer-Zarnowitz regressions

    y = forecastEvaluationData;
    X = staggeredCBO(:,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
    mzMat(ii,:,1) = coeff;
    mzMat_se(ii,:,1) = se;
    mzMat_p(ii,:,1) = 2*tcdf(-abs((coeff-[0;1])./se),floor(B));

    % Coibion-Gorodnichenko regressions

    if ii < (maxHorizon+1)

        y = cboForecastErrors(:,ii);
        X = staggeredCBO(:,ii)-staggeredCBO(:,ii+1);
        keepInd = ~isnan(y)&~isnan(X);
        y = y(keepInd);
        X = X(keepInd);
        TCur = numel(y);
        S = ceil(1.3*sqrt(TCur));
        B = TCur/S;

        [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');
        cgMat(ii,:,1) = coeff;
        cgMat_se(ii,:,1) = se;
        cgMat_p(ii,:,1) = 2*tcdf(-abs(coeff./se),floor(B));

    end

end

[biasMat(1:5)';biasMat_se(1:5)';biasMat_p(1:5)']
[arMat(1:5,2)';arMat_se(1:5,2)';arMat_p(1:5,2)']
[mzMat(1:5,2)';mzMat_se(1:5,2)';mzMat_p(1:5,2)']
[cgMat(1:5,2)';cgMat_se(1:5,2)';cgMat_p(1:5,2)']