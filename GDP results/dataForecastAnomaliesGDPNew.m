%% Load data and correct timing of forecasts

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

save('Gdp_Regression_Stats_New.mat','biasMat','biasMat_se','biasMat_p',...
    'arMat','arMat_se','arMat_p','mzMat','mzMat_se','mzMat_p',...
    'cgMat','cgMat_se','cgMat_p')