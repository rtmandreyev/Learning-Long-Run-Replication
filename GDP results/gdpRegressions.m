clear
clc

%% Load data

load('gdpResults_AllCoefficients_Levels5.mat')
load('Gdp_Data.mat')
load('Gdp_Regression_Stats.mat')

meanForecasts  = squeeze(mean(yGrowthForecastsADraws,3,'omitnan'));

meanForecasts(isnan(cboData)) = NaN;
staggeredMeanForecasts = meanForecasts;
maxHorizon = 10;
for ii = 0:maxHorizon

    staggeredMeanForecasts(:,ii+1) = circshift(meanForecasts(:,ii+1),ii);
    staggeredMeanForecasts(1:ii,ii+1) = NaN;

end

forecastEvaluationData(end) = NaN;
meanForecastErrors = forecastEvaluationData - staggeredMeanForecasts;

T = size(cboData,1);

biasMatModel = NaN(maxHorizon+1,1);
biasMatModel_se = biasMatModel;
biasMatModel_p = biasMatModel;

arMatModel = NaN(maxHorizon+1,2);
arMatModel_se = arMatModel;
arMatModel_p = arMatModel;

mzMatModel = NaN(maxHorizon+1,2);
mzMatModel_se = mzMatModel;
mzMatModel_p = mzMatModel;

cgMatModel = NaN(maxHorizon+1,2);
cgMatModel_se = cgMatModel;
cgMatModel_p = cgMatModel;

for ii = 1:(maxHorizon+1)

    % Bias

    X = ones(T-ii+1,1);
    y = meanForecastErrors(ii:end,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off','Intercept',false);

    biasMatModel(ii) = coeff;
    biasMatModel_se(ii) = se;
    biasMatModel_p(ii) = 2*tcdf(-abs(coeff./se),floor(B));

    % Autocorrelation of forecast errors

    X = meanForecastErrors(1:end-ii,ii);
    y = meanForecastErrors(ii+1:end,ii);
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');

    arMatModel(ii,:) = coeff;
    arMatModel_se(ii,:) = se;
    arMatModel_p(ii,:) = 2*tcdf(-abs(coeff./se),floor(B));

    % Mincer-Zarnowitz regressions

    X = staggeredMeanForecasts(:,ii);
    y = forecastEvaluationData;
    keepInd = ~isnan(y)&~isnan(X);
    y = y(keepInd);
    X = X(keepInd);
    TCur = numel(y);
    S = ceil(1.3*sqrt(TCur));
    B = TCur/S;

    [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');

    mzMatModel(ii,:) = coeff;
    mzMatModel_se(ii,:) = se;
    mzMatModel_p(ii,:) = 2*tcdf(-abs((coeff-[0;1])./se),floor(B));

    % Coibion-Gorodnichenko regressions

    if ii < (maxHorizon+1)

        X = staggeredMeanForecasts(:,ii)-...
            staggeredMeanForecasts(:,ii+1);
        y = meanForecastErrors(:,ii);
        keepInd = ~isnan(y)&~isnan(X);
        y = y(keepInd);
        X = X(keepInd);
        TCur = numel(y);
        S = ceil(1.3*sqrt(TCur));
        B = TCur/S;

        [~,se,coeff] = hac(X,y,'Bandwidth',S,'Display','off');

        cgMatModel(ii,:) = coeff;
        cgMatModel_se(ii,:) = se;
        cgMatModel_p(ii,:) = 2*tcdf(-abs(coeff./se),floor(B));

    end

end

% model results

% For each of the tables below, the first row contains the coefficient
% estimates, the second row contains the standard errors and the third row
% contains the p-values

[biasMatModel(1:5)';biasMatModel_se(1:5)';biasMatModel_p(1:5)'] % bias regression results
[arMatModel(1:5,2)';arMatModel_se(1:5,2)';arMatModel_p(1:5,2)'] % autocorrelated errors regression results
[mzMatModel(1:5,2)';mzMatModel_se(1:5,2)';mzMatModel_p(1:5,2)'] % Mincer-Zarnowitz regression results
[cgMatModel(1:5,2)';cgMatModel_se(1:5,2)';cgMatModel_p(1:5,2)'] % Coibion-Gorodnichenko regression results


