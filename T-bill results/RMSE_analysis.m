clear
clc
%% Computing RMSE for the Data
load('Tbl_Data.mat')
load('Tbl_Regression_Stats.mat')

% BC forecast errors are also computed

staggeredSPF = spfForecastsQ;
staggeredBC1 = bcForecasts1Q;
staggeredBC2 = bcForecasts2Q;
staggeredBC3 = bcForecasts3Q;
maxHorizonBC = 7;
maxHorizonSPF = 4;
numSeries = 4;

for ii = 0:maxHorizonBC
    
    if ii <= maxHorizonSPF
        staggeredSPF(:,ii+1) = circshift(spfForecastsQ(:,ii+1),ii);
        staggeredSPF(1:ii,ii+1) = NaN;
    end
    
    staggeredBC1(:,ii+1) = circshift(bcForecasts1Q(:,ii+1),ii);
    staggeredBC1(1:ii,ii+1) = NaN;
    
    staggeredBC2(:,ii+1) = circshift(bcForecasts2Q(:,ii+1),ii);
    staggeredBC2(1:ii,ii+1) = NaN;
    
    staggeredBC3(:,ii+1) = circshift(bcForecasts3Q(:,ii+1),ii);
    staggeredBC3(1:ii,ii+1) = NaN;
    
end

spfForecastErrors = tblData - staggeredSPF;
bc1ForecastErrors = tblData - staggeredBC1;
bc2ForecastErrors = tblData - staggeredBC2;
bc3ForecastErrors = tblData - staggeredBC3;
%% Getting rid of NaNs and Computing data RMSE
 

SPF_Nan_Row_sum=sum(isnan(spfForecastErrors),2);
SPF_Nan_Row_ind=find(SPF_Nan_Row_sum==0);
% Index of row from where spfForecastErrors do not have nans 
SPF_Notnan_ind=SPF_Nan_Row_ind(1);

spfForecastErrors_notnan=spfForecastErrors(SPF_Notnan_ind:end,1:end);

% RMSE of the data by quarter of forecast; Includes the Now cast


%% Getting rid of the nowcast and computing SPF forecast overall


SPF_RMSE_Quarter= sqrt(mean(spfForecastErrors_notnan.^2))

spfForecastErrors_notnan_4=spfForecastErrors_notnan (:,2:5)
SPF_RMSE_Full = sqrt(mean(spfForecastErrors_notnan_4(:).^2))

%% Baseline analysis



seed = 15284;

% thetLB = [0   ; 0  ; 0  ; 0         ];
% thetUB = [0.95; 0.1; 0.5; sqrt(1/12)];
% thetGuess = [0.8028; 0.0021; 0.0373; 0.0653]';
% options = optimoptions(@particleswarm,'Display','iter','OutputFcn',@gibbsUC_outfun_ps,...
%     'UseParallel',true,'InitialSwarmMatrix',thetGuess);
% thetOpt = particleswarm(@(thet) gibbsUC_Fun(thet,seed),4,thetLB',...
%     thetUB',options)

% Baseline Analysis:

thetOpt = [0.6; (0.12)^2; 0.09; 0.08];
rhoParams = thetOpt(1:2);
gamMoments = thetOpt(3:4);
if size(gamMoments,1) < size(gamMoments,2)
    gamMoments = gamMoments';
end
gamFun = @(x) [(x(1)-1)/(sum(x)-2); sqrt(prod(x)/((sum(x)^2)*(sum(x)+1)))];
options = optimoptions(@fminunc,'Display','off');
gamParams = fminunc(@(x) 1000*sum(((gamFun(x)-gamMoments)./gamMoments).^2),[2;2],options);
gamParams=[2.3,19.7];
gamParams_baseline=gamParams
%%
sig2Params = [1.25; 0.5625];
%%
burnIn = 50000;
B = 50000;
stepSize = 4;
%%
[rhoDraws,gamDraws,sigDraws,muRealTimeDraws,xRealTimeDraws,muSmoothedDraws,...
    xSmoothedDraws,yForecasts] = ...
    recursiveGibbsUC(tblData,rhoParams,gamParams,sig2Params,burnIn,B,...
    estimationStart,stepSize,(tblData>0.25),seed);

meanForecasts  = squeeze(mean(yForecasts,3,'omitnan'));
meanForecastsCopy = meanForecasts;

%% Compute model implied forecast errors

maxHorizon = 7;
meanForecastsCopy(isnan(spfForecastsQ(:,1)),:) = NaN;
staggeredModelForecasts = meanForecastsCopy(:,1:(maxHorizon+1));

for ii = 0:maxHorizon
    
    staggeredModelForecasts(:,ii+1) = circshift(meanForecastsCopy(:,ii+1),ii);
    staggeredModelForecasts(1:ii,ii+1) = NaN;
    
end

modelForecastErrors = tblData - staggeredModelForecasts;
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

for ii = 0:maxHorizon
    
    % Bias
    
    temp = regstats2(modelForecastErrors(ii+1:end,ii+1),ones(T-ii,1),...
        'onlydata',{'beta','hac'});
    biasMatModel(ii+1) = temp.beta;
    biasMatModel_se(ii+1) = temp.hac.se;
    biasMatModel_p(ii+1) = temp.hac.pval;
    
    % Autocorrelation of forecast errors
    
    temp = regstats2(modelForecastErrors(ii+1:end,ii+1),modelForecastErrors(1:end-ii,ii+1),...
        'linear',{'beta','hac'});
    if numel(temp.beta) == 1
        clear temp
        temp.beta = NaN(1,2);
        temp.hac.se = NaN(1,2);
        temp.hac.pval = NaN(1,2);
    end
    arMatModel(ii+1,:) = temp.beta;
    arMatModel_se(ii+1,:) = temp.hac.se;
    arMatModel_p(ii+1,:) = temp.hac.pval;
    
    % Mincer-Zarnowitz regressions
    
    temp = regstats2(tblData,staggeredModelForecasts(:,ii+1),'linear',{'beta','hac'});
    mzMatModel(ii+1,:) = temp.beta;
    mzMatModel_se(ii+1,:) = temp.hac.se;
    mzMatModel_p(ii+1,:) = [temp.hac.pval(1),2*normcdf(-abs((temp.beta(2)-1)/temp.hac.se(2)))];
    
    % Coibion-Gorodnichenko regressions
    
    if ii > 0
        
        temp = regstats2(modelForecastErrors(:,ii),staggeredModelForecasts(:,ii)-...
            staggeredModelForecasts(:,ii+1),'linear',{'beta','hac'});
        cgMatModel(ii,:) = temp.beta;
        cgMatModel_se(ii,:) = temp.hac.se;
        cgMatModel_p(ii,:) = temp.hac.pval;
        
    end
    
end

[biasMatModel, biasMatModel_se, biasMatModel_p]'
[arMatModel(:,2), arMatModel_se(:,2), arMatModel_p(:,2)]'
[mzMatModel(:,2), mzMatModel_se(:,2), mzMatModel_p(:,2)]'
[cgMatModel(:,2), cgMatModel_se(:,2), cgMatModel_p(:,2)]'


%% Computing RMSE for model forecast errors in the baseline case:

modelForecastErrors_SPFind_baseline = modelForecastErrors(SPF_Notnan_ind:end,2:5)

Baseline_RMSE_quarter=sqrt(mean(modelForecastErrors_SPFind_baseline.^2))
Baseline_RMSE_full=sqrt(mean(modelForecastErrors_SPFind_baseline(:).^2))

%% Defining the Theta that will be used to simulate the model with more dispersed priors

thet_dp1 = [0.6; (0.31)^2; 0.09; 0.10];




%% Running alternate scenario 1 and computing RMSE


rhoParams = thet_dp1(1:2);
gamMoments = thet_dp1(3:4);
if size(gamMoments,1) < size(gamMoments,2)
    gamMoments = gamMoments';
end
gamFun = @(x) [(x(1)-1)/(sum(x)-2); sqrt(prod(x)/((sum(x)^2)*(sum(x)+1)))];
options = optimoptions(@fminunc,'Display','off');
gamParams = fminunc(@(x) 1000*sum(((gamFun(x)-gamMoments)./gamMoments).^2),[2;2],options);
gamParams=[1.13,2.9];
gamParams_dp1=gamParams;
%%
sig2Params = [1.25; 0.5625];
burnIn = 50000;
B = 50000;
stepSize = 4;
[rhoDraws,gamDraws,sigDraws,muRealTimeDraws,xRealTimeDraws,muSmoothedDraws,...
    xSmoothedDraws,yForecasts] = ...
    recursiveGibbsUC_sh(tblData,rhoParams,gamParams,sig2Params,burnIn,B,...
    estimationStart,stepSize,(tblData>0.25),seed);

meanForecasts  = squeeze(mean(yForecasts,3,'omitnan'));
meanForecastsCopy = meanForecasts;


% Computing forecast errors
maxHorizon = 7;
meanForecastsCopy(isnan(spfForecastsQ(:,1)),:) = NaN;
staggeredModelForecasts = meanForecastsCopy(:,1:(maxHorizon+1));

for ii = 0:maxHorizon
    
    staggeredModelForecasts(:,ii+1) = circshift(meanForecastsCopy(:,ii+1),ii);
    staggeredModelForecasts(1:ii,ii+1) = NaN;
    
end

modelForecastErrors = tblData - staggeredModelForecasts;

% Computing the forecast errors for the appropriate periods in alternate
% scenario 1

modelForecastErrors_SPFind_DP1 = modelForecastErrors(SPF_Notnan_ind:end,2:5);
Dispersed_P1_RMSE_quarter=sqrt(mean(modelForecastErrors_SPFind_DP1.^2))
Dispersed_P1_RMSE_full=sqrt(mean(modelForecastErrors_SPFind_DP1(:).^2))


%%
Baseline_RMSE_quarter=sqrt(mean(modelForecastErrors_SPFind_baseline.^2))
Baseline_RMSE_full=sqrt(mean(modelForecastErrors_SPFind_baseline(:).^2))


%% Table F.5
%Baseline RMSE ratio
[Baseline_RMSE_quarter./(SPF_RMSE_Quarter(2:5)) , Baseline_RMSE_full/SPF_RMSE_Full]


%dispersed priors RMSE ratio
[Dispersed_P1_RMSE_quarter./SPF_RMSE_Quarter(2:5) , Dispersed_P1_RMSE_full/SPF_RMSE_Full]