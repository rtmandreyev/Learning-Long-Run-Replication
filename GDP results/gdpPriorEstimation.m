clear
clc

load('Gdp_Data.mat')
load('Gdp_Regression_Stats.mat')

burnIn = 25000;
B = 50000;
stepSize = 4;
seed = 15284;

% thetGuess = [0.5991;0.1315^2;-0.8344;0.0122^2;10.6146; 9.1436];
% thetGuess = [0.8259; 1e-4; -0.8048; 1e-4; 44.1410; 100];
%thetGuess = [0.5211; 0.01; -0.9982; 0.0101; 15.9760; 42.0606];
% thetGuess = [0.6961; 0.01; -0.9902; 0.0102; 23.7525; 18.2262];
% options = optimoptions(@particleswarm,'Display','iter','OutputFcn',@gibbsUC_GDP_outfun_ps,...
%     'UseParallel',true,'InitialSwarmMatrix',thetGuess');
% thetLB = [0.5; 1e-4; -1; 1e-4;   1;   1];
% thetUB = [  1; 0.25;  0; 0.25; 100; 100];
% 
% thetOpt = particleswarm(@(thet) gibbsUC_GDP_Fun(thet,seed),6,thetLB',...
%     thetUB',options);

%thetOpt = [0.5343; 0.012243; -0.6830; 1.4892e-4; 22.0006; 50.4930];
%thetOpt = [0.8724; 1e-4; -0.8757; 1e-4; 1.9115; 3.5255];
%thetOpt = [0.5991;0.1315^2;-0.8344;0.0122^2;10.6146; 9.1436];
%thetOpt = [0.6703; 0.01; -0.9982; 0.01; 34.7091; 28.0402];
thetOpt = [0.8888; 0.0028; -0.9129; 0.0025; 45.9721; 39.3704];
rhoParams = thetOpt(1:2);
rho2Params = thetOpt(3:4);
gamParams = thetOpt(5:6);
muParams = [0.01; 0.01^2];
%muParams = thetOpt(7:8);

%muParams = [0.0125; 0.015^2];
sig2Fun = @(x) [x(2)/(x(1)+1); x(2)/sqrt((x(1)-1)^2*(x(1)-2))];
sig2Moments = [0.015^2; 0.01^2];
sig2Params = [sig2Moments(1)^2/(sig2Moments(2)^2)+2; NaN];
sig2Params(2) = sig2Moments(1)*(sig2Params(1)-1);

[muDraws,gamDraws,rho1Draws,rho2Draws,sigDraws,zRealTimeDraws,xRealTimeDraws,...
    zSmoothedDraws,xSmoothedDraws,~,yGrowthForecastsADraws] = ...
    recursiveGibbsUC_GDP(log(gdpData),muParams,gamParams,rhoParams,...
    rho2Params,sig2Params,burnIn,B,estimationStart,vintageStart,...
    stepSize,true(size(gdpData,1),1),seed,1);

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
    
    temp = regstats2(meanForecastErrors(ii:end,ii),ones(T-ii+1,1),...
        'onlydata',{'beta','hac'});
    biasMatModel(ii) = temp.beta;
    biasMatModel_se(ii) = temp.hac.se;
    biasMatModel_p(ii) = temp.hac.pval;
    
    % Autocorrelation of forecast errors
    
    temp = regstats2(meanForecastErrors(ii+1:end,ii),meanForecastErrors(1:end-ii,ii),...
        'linear',{'beta','hac'});
    if numel(temp.beta) == 1
        temp.beta = NaN(1,2);
        temp.hac.se = NaN(1,2);
        temp.hac.pval = NaN(1,2);
    end
    arMatModel(ii,:) = temp.beta;
    arMatModel_se(ii,:) = temp.hac.se;
    arMatModel_p(ii,:) = temp.hac.pval;
    
    % Mincer-Zarnowitz regressions
    
    temp = regstats2(forecastEvaluationData,staggeredMeanForecasts(:,ii),'linear',{'beta','hac'});
    mzMatModel(ii,:) = temp.beta;
    mzMatModel_se(ii,:) = temp.hac.se;
    mzMatModel_p(ii,:) = [temp.hac.pval(1),2*normcdf(-abs((temp.beta(2)-1)/temp.hac.se(2)))];
    
    % Coibion-Gorodnichenko regressions
    
    if ii < (maxHorizon+1)
        
        temp = regstats2(meanForecastErrors(:,ii),staggeredMeanForecasts(:,ii)-...
            staggeredMeanForecasts(:,ii+1),'linear',{'beta','hac'});
        cgMatModel(ii,:) = temp.beta;
        cgMatModel_se(ii,:) = temp.hac.se;
        cgMatModel_p(ii,:) = temp.hac.pval;
        
    end
    
end

evalInd = 1:5;
tempBias = biasMat(evalInd) - biasMatModel(evalInd);
tempBias_se = biasMat_se(evalInd) - biasMatModel_se(evalInd);
biasDenom = biasMatModel(evalInd);
tempAR = arMat(evalInd,:) - arMatModel(evalInd,:);
tempAR_se = arMat_se(evalInd,:) - arMatModel_se(evalInd,:);
arDenom = arMat(evalInd,:);
tempMZ = mzMat(evalInd,:) - mzMatModel(evalInd,:);
tempMZ_se = mzMat_se(evalInd,:) - mzMatModel_se(evalInd,:);
mzDenom = mzMat(evalInd,:);
tempCG = cgMat(evalInd,:) - cgMatModel(evalInd,:);
tempCG_se = cgMat_se(evalInd,:) - cgMatModel_se(evalInd,:);
cgDenom = cgMat(evalInd,:);

% obj = sum([tempBias./biasDenom; tempAR(:)./arDenom(:);...
%     tempMZ(:)./mzDenom(:); tempCG(:)./cgDenom(:)].^2,'omitnan')

%obj = sum([tempBias; tempAR(:,2); tempMZ(:,2); tempCG(:,2)].^2,'omitnan')

obj = sum([tempBias; tempAR(:); tempMZ(:); tempCG(:)].^2,'omitnan')
%obj = sum([tempBias; tempAR(:,2); tempMZ(:,2); tempCG(:,2); ...
%        tempBias_se; tempAR_se(:,2); tempMZ_se(:,2); tempCG_se(:,2)].^2,'omitnan')

% obj = [tempBias; tempAR(:,2); tempMZ(:,2); tempCG(:,2)]'*...
%         (diag(1./([tempBias_se; tempAR_se(:,2); tempMZ_se(:,2); tempCG_se(:,2)].^2)))*...
%         [tempBias; tempAR(:,2); tempMZ(:,2); tempCG(:,2)]

save('gdpResults_AllCoefficients_Levels5.mat','biasMatModel',...
    'biasMatModel_se','biasMatModel_p', ...
    'arMatModel','arMatModel_se',...
    'arMatModel_p','mzMatModel','mzMatModel_se','mzMatModel_p',...
    'cgMatModel','cgMatModel_se','cgMatModel_p','rho1Draws','rho2Draws',...
    'gamDraws','muDraws','sigDraws','zRealTimeDraws','xRealTimeDraws',...
    'zSmoothedDraws','xSmoothedDraws','yGrowthForecastsADraws','rhoParams',...
    'rho2Params','gamParams','muParams','sig2Params')