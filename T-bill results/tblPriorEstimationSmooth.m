clear
clc

load('Tbl_Data.mat')
load('Tbl_Regression_Stats.mat')
%%
seed = 15284;

thetLB = [0   ; 0  ; 0  ; 0         ];
thetUB = [0.95; 0.1; 0.5; sqrt(1/12)];
% thetGuess = [0.6535; 0.0105; 0.0539; 0.0788]';
% options = optimoptions(@particleswarm,'Display','iter','OutputFcn',@gibbsUC_outfun_ps,...
%    'UseParallel',true,'InitialSwarmMatrix',thetGuess,'HybridFcn',@fmincon);
% thetOpt = particleswarm(@(thet) gibbsUC_Fun(thet,seed),4,thetLB',...
%    thetUB',options)
% thetOpt = fmincon
% (@(thet) gibbsUC_Fun(thet,seed),thetOpt,[],[],[],[],...
%     thetLB,thetUB,[],optimoptions(@fmincon,'Display','iter',...
%     'UseParallel',true))


thetOpt = [0.7605957; 0.0048179; 0.0117594; 0.0809342];
% parameters above are: rho_mean = 0.76, rho_variance=0.0048179

%thetOpt = [0.6037; 0.0144; 0.0654; 0.0638];
%thetOpt = [0.3397; 0.0233; 0.0301; 0.0532];
rhoParams = thetOpt(1:2);
gamMoments = thetOpt(3:4);
if size(gamMoments,1) < size(gamMoments,2)
    gamMoments = gamMoments';
end
gamFun = @(x) [(x(1)-1)/(sum(x)-2); sqrt(prod(x)/((sum(x)^2)*(sum(x)+1)))];
options = optimoptions(@fminunc,'Display','off');
gamParams = fminunc(@(x) 1000*sum(((gamFun(x)-gamMoments)./gamMoments).^2),[2;2],options);

sig2Params = [1.25; 0.5625];
burnIn = 50000;
B = 75000;
stepSize = 1;
%%
[rhoDraws,gamDraws,sigDraws,muRealTimeDraws,xRealTimeDraws,muSmoothedDraws,...
    xSmoothedDraws,yForecasts] = ...
    recursiveGibbsUC(tblData,rhoParams,gamParams,sig2Params,burnIn,B,...
    estimationStart,stepSize,(tblData>0.25),seed);

meanForecasts  = squeeze(mean(yForecasts,3,'omitnan'));
meanForecastsCopy = meanForecasts;

%% Compute model implied forecast anomalies

maxHorizon = 7;
meanForecastsCopy(isnan(spfForecastsQ(:,1)),:) = NaN;
staggeredModelForecasts = meanForecastsCopy(:,1:(maxHorizon+1));

for ii = 0:maxHorizon
    
    staggeredModelForecasts(:,ii+1) = circshift(meanForecastsCopy(:,ii+1),ii);
    staggeredModelForecasts(1:ii,ii+1) = NaN;
    
end

modelForecastErrors = tblData - staggeredModelForecasts;

%% Compute forecast anomalies

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

%% Compute forecast anomalies across sub-samples

sample1Ind = (1:155)';
sample2Ind = (156:T)';

% Compute forecast anomalies

biasMatModelSplit = NaN(maxHorizon+1,2);
biasMatModelSplit_se = biasMatModelSplit;
biasMatModelSplit_p = biasMatModelSplit;

arMatModelSplit = NaN(maxHorizon+1,2,2);
arMatModelSplit_se = arMatModelSplit;
arMatModelSplit_p = arMatModelSplit;

mzMatModelSplit = NaN(maxHorizon+1,2,2);
mzMatModelSplit_se = mzMatModelSplit;
mzMatModelSplit_p = mzMatModelSplit;

cgMatModelSplit = NaN(maxHorizon+1,2,2);
cgMatModelSplit_se = cgMatModelSplit;
cgMatModelSplit_p = cgMatModelSplit;

for split = 1:2
    
    if split == 1
        sampleInd = sample1Ind;
    else
        sampleInd = sample2Ind;
    end
    
    for ii = 0:maxHorizon
        
        % Bias
        
        temp = regstats2(modelForecastErrors(sampleInd(ii+1:end),ii+1),...
            ones(numel(sampleInd)-ii,1),'onlydata',{'beta','hac'});
        biasMatModelSplit(ii+1,split) = temp.beta;
        biasMatModelSplit_se(ii+1,split) = temp.hac.se;
        biasMatModelSplit_p(ii+1,split) = temp.hac.pval;
        
        % Autocorrelation of forecast errors
        
        temp = regstats2(modelForecastErrors(sampleInd(ii+1:end),ii+1),...
            modelForecastErrors(sampleInd(1:end-ii),ii+1),'linear',{'beta','hac'});
        if numel(temp.beta) == 1
            clear temp
            temp.beta = NaN(1,2);
            temp.hac.se = NaN(1,2);
            temp.hac.pval = NaN(1,2);
        end
        arMatModelSplit(ii+1,:,split) = temp.beta;
        arMatModelSplit_se(ii+1,:,split) = temp.hac.se;
        arMatModelSplit_p(ii+1,:,split) = temp.hac.pval;
        
        % Mincer-Zarnowitz regressions
        
        temp = regstats2(tblData(sampleInd),staggeredModelForecasts(sampleInd,ii+1),...
            'linear',{'beta','hac'});
        mzMatModelSplit(ii+1,:,split) = temp.beta;
        mzMatModelSplit_se(ii+1,:,split) = temp.hac.se;
        mzMatModelSplit_p(ii+1,:,split) = [temp.hac.pval(1),...
            2*normcdf(-abs((temp.beta(2)-1)/temp.hac.se(2)))];
        
        % Coibion-Gorodnichenko regressions
        
        if ii > 0
            
            temp = regstats2(modelForecastErrors(sampleInd,ii),...
                staggeredModelForecasts(sampleInd,ii)-...
                staggeredModelForecasts(sampleInd,ii+1),'linear',...
                {'beta','hac'});
            cgMatModelSplit(ii,:,split) = temp.beta;
            cgMatModelSplit_se(ii,:,split) = temp.hac.se;
            cgMatModelSplit_p(ii,:,split) = temp.hac.pval;
            
        end
        
    end
    
    [biasMatModelSplit(:,split), biasMatModelSplit_se(:,split), biasMatModelSplit_p(:,split)]'
    [squeeze(arMatModelSplit(:,2,split)), squeeze(arMatModelSplit_se(:,2,split)), squeeze(arMatModelSplit_p(:,2,split))]'
    [squeeze(mzMatModelSplit(:,2,split)), squeeze(mzMatModelSplit_se(:,2,split)), squeeze(mzMatModelSplit_p(:,2,split))]'
    [squeeze(cgMatModelSplit(:,2,split)), squeeze(cgMatModelSplit_se(:,2,split)), squeeze(cgMatModelSplit_p(:,2,split))]'
    
end

%% Compute expectations hypothesis violations

maxHorizon = 40;
breakInd = 156;

modelYields = cumsum(meanForecasts,2)./(1:maxHorizon);

mVec = 1; % short horizon
numM = numel(mVec);
nVec = [2 3 4 8 12 20 40]; % long horizon
numN = numel(nVec);
csConventionalMatModel = NaN(numM,numN);
csConventionalMatModel_se = csConventionalMatModel;
csConventionalMatModel_p = csConventionalMatModel;
csContrarianMatModel = NaN(numM,numN);
csContrarianMatModel_se = csContrarianMatModel;
csContrarianMatModel_p = csContrarianMatModel;

for mInd = 1:numM
    
    m = mVec(mInd);
    
    for nInd = 1:numN
        
        n = nVec(nInd);
        
        yieldInd = (sum(isnan(yieldData(:,n)))+1:T);
        %yieldInd = yieldInd(yieldInd>=breakInd);
        %yieldInd = yieldInd(yieldInd<breakInd);
        
        if ~mod(n,m)
            
            k = n/m;
            
            TReg1 = numel(yieldInd) - (k-1)*m;
            xReg = modelYields(yieldInd(1:TReg1),n) - modelYields(yieldInd(1:TReg1),m);
            yReg1 = NaN(size(xReg));
            for t = 1:TReg1
                
                yReg1(t) = mean(modelYields(yieldInd(t:m:t+m*(k-1)),m)) - ...
                    modelYields(yieldInd(t),m);
                
            end
            
            temp = regstats2(yReg1,xReg,'linear',{'beta','hac'});
            csConventionalMatModel(mInd,nInd) = temp.beta(2);
            csConventionalMatModel_se(mInd,nInd) = temp.hac.se(2);
            csConventionalMatModel_p(mInd,nInd) = 2*normcdf(-abs((temp.beta(2)-1)./temp.hac.se(2)));
            
            TReg2 = numel(yieldInd) - m;
            xReg = (m/(n-m))*(modelYields(yieldInd(1:TReg2),n) - modelYields(yieldInd(1:TReg2),m));
            yReg2 = modelYields(yieldInd(1+m:end),n-m) - modelYields(yieldInd(1:TReg2),n);
            
            temp = regstats2(yReg2,xReg,'linear',{'beta','hac'});
            csContrarianMatModel(mInd,nInd) = temp.beta(2);
            csContrarianMatModel_se(mInd,nInd) = temp.hac.se(2);
            csContrarianMatModel_p(mInd,nInd) = 2*normcdf(-abs((temp.beta(2)-1)./temp.hac.se(2)));
            
        end
        
    end
end

[csConventionalMatModel; csConventionalMatModel_se; csConventionalMatModel_p]
[csContrarianMatModel; csContrarianMatModel_se; csContrarianMatModel_p]

%% Compute expectations hypothesis violations across sub-samples

csConventionalMatModelSplit = NaN(numM,numN,2);
csConventionalMatModelSplit_se = csConventionalMatModelSplit;
csConventionalMatModelSplit_p = csConventionalMatModelSplit;
csContrarianMatModelSplit = NaN(numM,numN,2);
csContrarianMatModelSplit_se = csContrarianMatModelSplit;
csContrarianMatModelSplit_p = csContrarianMatModelSplit;

for split = 1:2
    
    for mInd = 1:numM
        
        m = mVec(mInd);
        
        for nInd = 1:numN
            
            n = nVec(nInd);
            
            yieldInd = (sum(isnan(yieldData(:,n)))+1:T);
            if split == 1
                yieldInd = yieldInd(yieldInd<breakInd);
            else
                yieldInd = yieldInd(yieldInd>=breakInd);
            end
            
            if ~mod(n,m)
                
                k = n/m;
                
                TReg1 = numel(yieldInd) - (k-1)*m;
                xReg = modelYields(yieldInd(1:TReg1),n) - modelYields(yieldInd(1:TReg1),m);
                yReg1 = NaN(size(xReg));
                for t = 1:TReg1
                    
                    yReg1(t) = mean(modelYields(yieldInd(t:m:t+m*(k-1)),m)) - ...
                        modelYields(yieldInd(t),m);
                    
                end
                
                temp = regstats2(yReg1,xReg,'linear',{'beta','hac'});
                csConventionalMatModelSplit(mInd,nInd,split) = temp.beta(2);
                csConventionalMatModelSplit_se(mInd,nInd,split) = temp.hac.se(2);
                csConventionalMatModelSplit_p(mInd,nInd,split) = ...
                    2*normcdf(-abs((temp.beta(2)-1)./temp.hac.se(2)));
                
                TReg2 = numel(yieldInd) - m;
                xReg = (m/(n-m))*(modelYields(yieldInd(1:TReg2),n) - ...
                    modelYields(yieldInd(1:TReg2),m));
                yReg2 = modelYields(yieldInd(1+m:end),n-m) - ...
                    modelYields(yieldInd(1:TReg2),n);
                
                temp = regstats2(yReg2,xReg,'linear',{'beta','hac'});
                csContrarianMatModelSplit(mInd,nInd,split) = temp.beta(2);
                csContrarianMatModelSplit_se(mInd,nInd,split) = temp.hac.se(2);
                csContrarianMatModelSplit_p(mInd,nInd,split) = ...
                    2*normcdf(-abs((temp.beta(2)-1)./temp.hac.se(2)));
                
            end
            
        end
    end
    
    [squeeze(csConventionalMatModelSplit(:,:,split)); ...
        squeeze(csConventionalMatModelSplit_se(:,:,split)); ...
        squeeze(csConventionalMatModelSplit_p(:,:,split))]
    [squeeze(csContrarianMatModelSplit(:,:,split)); ...
        squeeze(csContrarianMatModelSplit_se(:,:,split)); ...
        squeeze(csContrarianMatModelSplit_p(:,:,split))]
    
end

% obj = sum([biasMat(2:end,1)-biasMatModel(2:end);biasMat_se(2:end,1)-...
%     biasMatModel_se(2:end);arMat(2:end,2,1)-arMatModel(2:end,2);...
%     arMat_se(2:end,2,1)-arMatModel_se(2:end,2);mzMat(2:end,2,1)-...
%     mzMatModel(2:end,2);mzMat_se(2:end,2,1)-mzMatModel_se(2:end,2);...
%     cgMat(2:end,2,1)-cgMatModel(2:end,2);cgMat_se(2:end,2,1)-...
%     cgMatModel_se(2:end,2); csConventionalMat'-csConventionalMatModel';...
%     csConventionalMat_se'-csConventionalMatModel_se';csContrarianMat'-...
%     csContrarianMatModel';csContrarianMat_se'-csContrarianMatModel_se'].^2,...
%     'omitnan')

obj = sum([biasMat(2:end,1)-biasMatModel(2:end);arMat(2:end,2,1)-arMatModel(2:end,2);...
    mzMat(2:end,2,1)-mzMatModel(2:end,2);...
    cgMat(2:end,2,1)-cgMatModel(2:end,2); csConventionalMat'-csConventionalMatModel';...
    csContrarianMat'-csContrarianMatModel'].^2,...
    'omitnan')

yieldModel = modelYields;
yForecastsShort = yForecasts(:,1:10,:);

save('tblResults_SlopeCoefficients_Level_SmoothL.mat','biasMatModel',...
    'biasMatModel_se','biasMatModel_p','arMatModel','arMatModel_se',...
    'arMatModel_p','mzMatModel','mzMatModel_se','mzMatModel_p',...
    'cgMatModel','cgMatModel_se','cgMatModel_p','rhoDraws','gamDraws',...
    'sigDraws','muRealTimeDraws','xRealTimeDraws','muSmoothedDraws',...
    'xSmoothedDraws','yForecastsShort','csConventionalMatModel',...
    'csConventionalMatModel_se','csConventionalMatModel_p',...
    'csContrarianMatModel','csContrarianMatModel_se',...
    'csContrarianMatModel_p','yieldModel','rhoParams','gamParams',...
    'sig2Params','biasMatModelSplit','biasMatModelSplit_se',...
    'biasMatModelSplit_p','arMatModelSplit','arMatModelSplit_se',...
    'arMatModelSplit_p','mzMatModelSplit','mzMatModelSplit_se',...
    'mzMatModelSplit_p','cgMatModelSplit','cgMatModelSplit_se',...
    'cgMatModelSplit_p','csConventionalMatModelSplit','csConventionalMatModelSplit_se',...
    'csConventionalMatModelSplit_p','csContrarianMatModelSplit',...
    'csContrarianMatModelSplit_se','csContrarianMatModelSplit_p')