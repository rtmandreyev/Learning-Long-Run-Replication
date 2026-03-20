
clear,clc;
data = 'New_DownwardBiasedPriors';
savename = strcat('simulationGibbs_',data,'.mat');


%% True model parameters

rho = 0.95;
gam = 0.3;
sig = 0.5;
trendInd = 0;

%% Prior parameters (modify to match existing results)

rhoParams = [0.95; 0.1^2]; % mean and variance of Normal prior for rho
gamMoments = [0.01; 0.01]; % mode and standard deviation of Beta prior for gamma
gamFun = @(x) [(x(1)-1)/(sum(x)-2); sqrt(prod(x)/((sum(x)^2)*(sum(x)+1)))]; % moments as function of Beta distribution parameters
gamParams = fminunc(@(x) 1000*sum((gamFun(x)-gamMoments).^2),[2;2]); % solve for Beta distribution parameters associated with moments
gamParams = [9.052;19.788];
sig2Params = [1.25; 0.5625]; % Inverse Gamma distribution prior parameters for conditional variance, do not change

invGamPdf = @(x,a,b) ((b^a)/gamma(a))*(x.^(-a-1)).*exp(-b./x).*(x>0);

% Distributions
rhoPrior = @(x) normpdf(x,rhoParams(1),sqrt(rhoParams(2)));
gamPrior = @(x) betapdf(x,gamParams(1),gamParams(2));
sigPrior = @(x) invGamPdf(x,sig2Params(1),sig2Params(2));

priorPlot1 = linspace(0,1,1000)';
priorPlot2 = linspace(0,2,1000)';
figure
subplot(3,1,1)
plot(priorPlot1,rhoPrior(priorPlot1),'k-','LineWidth',3)
xline(rho,'Color','#808080','LineWidth',3);
xt = xticks;
title('\rho')
set(gca,'FontSize',20)
axis tight

subplot(3,1,2)
plot(priorPlot1,gamPrior(priorPlot1),'k-','LineWidth',3)
xline(gam,'Color','#808080','LineWidth',3);
title('\gamma')
set(gca,'FontSize',20)
axis tight

subplot(3,1,3)
plot(priorPlot2,sigPrior(priorPlot2),'k-','LineWidth',3)
xline(sig^2,'Color','#808080','LineWidth',3);
title('\sigma^2')
set(gca,'FontSize',20)
axis tight

figurename =  strcat('Priors/simulationGibbs_',data);
%saveas(gcf,figurename,'epsc')



%% Analysis of Data

clear,clc;
load('Tbl_Data.mat','tblData','estimationStart','T','spfForecastsQ','yieldData')
load('Tbl_Stats_Final.mat')
data = 'New_DownwardBiasedPriors';
datafile = strcat('simulationGibbs_',data,'.mat');
load(datafile);
tblSamples = muSamples + xSamples;
circFlag = 0;
estimationStart = 42; %dateVec(42) = 1961Q3, start of estimation

%% Plot Parameters

%True Params
thetaTrue = NaN(3,T);
thetaTrue(1,:) = rho*ones(T,1);
thetaTrue(2,:) = gam*ones(T,1);
thetaTrue(3,:) = sig*ones(T,1);

%Simulated Params
thetaSim = NaN(3,T,B);
thetaSim(1,:,:) = rhoMeanSim;
thetaSim(2,:,:) = gamMeanSim;
thetaSim(3,:,:) = sigMeanSim;

maxHorizon = 7;

% plot parameters, with 95% CI

close all;
figure('units','normalized','outerposition',[0 0 1 1])

dateVec = datetime(1951,1*3 + 1,1) + calquarters(0:T-1);

subplot(3,1,1)
hold on;
plot(dateVec,thetaTrue(1,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(1,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(1,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(1,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\rho')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

subplot(3,1,2)
hold on;
plot(dateVec,thetaTrue(2,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(2,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(2,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(2,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\gamma')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

subplot(3,1,3)
hold on;
plot(dateVec,thetaTrue(3,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(3,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(3,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(3,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\sigma_{\nu}')
axis tight
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

%figurename =  strcat('Params/simulationGibbs_',data);
%saveas(gcf,figurename,'epsc')






%% Compute regression statistics

% Structures of forecast rationality regressions

biasSim = NaN(maxHorizon+1,B);
biasSim_se = biasSim;
biasSim_p = biasSim;

arSim = NaN(maxHorizon+1,2,B);
arSim_se = arSim;
arSim_p = arSim;

mzSim = NaN(maxHorizon+1,2,B);
mzSim_se = mzSim;
mzSim_p = mzSim;

cgSim = NaN(maxHorizon+1,2,B);
cgSim_se = cgSim;
cgSim_p = cgSim;

% Structures for Campbell-Shiller regressions

mVec = 1; % short horizon
numM = numel(mVec);
nVec = [2 3 4 8 12 20 40]; % long horizon
numN = numel(nVec);

csConventionalSim = NaN(numM,numN,B);
csConventionalSim_se = csConventionalSim;
csConventionalSim_p = csConventionalSim;

csContrarianSim = NaN(numM,numN,B);
csContrarianSim_se = csContrarianSim;
csContrarianSim_p = csContrarianSim;

spfForecastsQSlice = spfForecastsQ(:,1);
for b = 1:B
    
    meanForecastsCopy = meanForecastsSim(:,:,b);
    meanForecastsCopy(isnan(spfForecastsQSlice),:) = NaN;
    staggeredModelForecasts = meanForecastsCopy(:,1:(maxHorizon+1));
    
    for ii = 0:maxHorizon
        
        staggeredModelForecasts(:,ii+1) = circshift(meanForecastsCopy(:,ii+1),ii);
        staggeredModelForecasts(1:ii,ii+1) = NaN;
        
    end
    
    
    modelForecastErrors = tblSamples(:,b) - staggeredModelForecasts;
    if all(isnan(modelForecastErrors))
        %Do Nothing
        modelForecastErrors;
    else
    
    
        % Forecasting anomaly regressions

        for ii = 0:maxHorizon

            % Bias

            temp = regstats2(modelForecastErrors(ii+1:end,ii+1),ones(T-ii,1),...
                'onlydata',{'beta','hac'});
            biasSim(ii+1,b) = temp.beta;
            try
            biasSim_se(ii+1,b) = temp.hac.se;
            catch
                temp.hac.se;
            end
            biasSim_p(ii+1,b) = temp.hac.pval;

            % Autocorrelation of forecast errors

            temp = regstats2(modelForecastErrors(ii+1:end,ii+1),modelForecastErrors(1:end-ii,ii+1),...
                'linear',{'beta','hac'});
            if numel(temp.beta) == 1
                clear temp
                temp.beta = NaN(1,2);
                temp.hac.se = NaN(1,2);
                temp.hac.pval = NaN(1,2);
            end
            arSim(ii+1,:,b) = temp.beta;
            arSim_se(ii+1,:,b) = temp.hac.se;
            arSim_p(ii+1,:,b) = temp.hac.pval;

            % Mincer-Zarnowitz regressions

            temp = regstats2(tblSamples(:,b),staggeredModelForecasts(:,ii+1),'linear',{'beta','hac'});
            mzSim(ii+1,:,b) = temp.beta;
            mzSim_se(ii+1,:,b) = temp.hac.se;
            mzSim_p(ii+1,:,b) = [temp.hac.pval(1),2*normcdf(-abs((temp.beta(2)-1)/temp.hac.se(2)))];

            % Coibion-Gorodnichenko regressions

            if ii > 0

                temp = regstats2(modelForecastErrors(:,ii),staggeredModelForecasts(:,ii)-...
                    staggeredModelForecasts(:,ii+1),'linear',{'beta','hac'});
                cgSim(ii,:,b) = temp.beta;
                cgSim_se(ii,:,b) = temp.hac.se;
                cgSim_p(ii,:,b) = temp.hac.pval;

            end

        end

        % Campbell-Shiller regressions

        modelYields = cumsum(meanForecastsSim,2)./(1:40);

        for mInd = 1:numM

            m = mVec(mInd);

            for nInd = 1:numN

                n = nVec(nInd);

                yieldInd = (sum(isnan(yieldData(:,n)))+1:T);

                if ~mod(n,m)

                    k = n/m;

                    TReg1 = numel(yieldInd) - (k-1)*m;
                    xReg = modelYields(yieldInd(1:TReg1),n,b) - modelYields(yieldInd(1:TReg1),m,b);
                    yReg1 = NaN(size(xReg));
                    for t = 1:TReg1

                        yReg1(t) = mean(modelYields(yieldInd(t:m:t+m*(k-1)),m,b)) - ...
                            modelYields(yieldInd(t),m,b);

                    end

                    temp = regstats2(yReg1,xReg,'linear',{'beta','hac'});
                    csConventionalSim(mInd,nInd,b) = temp.beta(2);
                    csConventionalSim_se(mInd,nInd,b) = temp.hac.se(2);
                    csConventionalSim_p(mInd,nInd,b) = 2*normcdf(-abs((temp.beta(2)-1)./temp.hac.se(2)));

                    TReg2 = numel(yieldInd) - m;
                    xReg = (m/(n-m))*(modelYields(yieldInd(1:TReg2),n,b) - modelYields(yieldInd(1:TReg2),m,b));
                    yReg2 = modelYields(yieldInd(1+m:end),n-m,b) - modelYields(yieldInd(1:TReg2),n,b);

                    temp = regstats2(yReg2,xReg,'linear',{'beta','hac'});
                    csContrarianSim(mInd,nInd,b) = temp.beta(2);
                    csContrarianSim_se(mInd,nInd,b) = temp.hac.se(2);
                    csContrarianSim_p(mInd,nInd,b) = 2*normcdf(-abs((temp.beta(2)-1)./temp.hac.se(2)));

                end     
            end
        end
    end 

    
end



%% Table 6 and 7 Regression Results 
% Autocorrelation results
[nanmean(reshape(arSim(2:5,2,:),[],B)'); nanstd(reshape(arSim(2:5,2,:),[],B)'); ECDF(reshape(arSim(2:5,2,:),[],B)',dataStats.ar(2:end,2)')]

% Mincer Zarnowitz results
[nanmean(reshape(mzSim(2:5,2,:),[],B)'); nanstd(reshape(mzSim(2:5,2,:),[],B)'); ECDF(reshape(mzSim(2:5,2,:),[],B)',dataStats.mz(2:end,2)')]

% Coibion Gorodnichenko results
[nanmean(reshape(cgSim(1:4,2,:),[],B)'); nanstd(reshape(cgSim(1:4,2,:),[],B)'); ECDF(reshape(cgSim(1:4,2,:),[],B)',dataStats.cg(1:4,2)')]

% Campbell Shiller Conventional results

[nanmean(permute(csConventionalSim,[3,2,1])); nanstd(permute(csConventionalSim,[3,2,1])); ECDF(permute(csConventionalSim,[3,2,1]),dataStats.csConventional')]

% Campbell Shiller Contrarian results

[nanmean(permute(csContrarianSim,[3,2,1]));nanstd(permute(csContrarianSim,[3,2,1])); ECDF(permute(csContrarianSim,[3,2,1]),dataStats.csContrarian')]



%% Figure 14 (when the downward biased runs data file is loaded in)
close all;
figure('units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)
hold on;
plot(dateVec,thetaTrue(1,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(1,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(1,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(1,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\rho')
axis square
set(gca,'FontSize',20)
ylim([0 1]);
xlim(dateVec([estimationStart T]));

subplot(1,2,2)
hold on;
plot(dateVec,thetaTrue(2,:),'Color','#808080','LineWidth',2);
plot(dateVec,nanmean(thetaSim(2,:,:),3)','k','LineWidth',2)
plot(dateVec,prctile(thetaSim(2,:,:),5,3)','k--','LineWidth',2)
plot(dateVec,prctile(thetaSim(2,:,:),95,3)','k--','LineWidth',2)
hold off;
title('\gamma')
axis square
set(gca,'FontSize',20)
ylim([0 .5]);
xlim(dateVec([estimationStart T]));

%figurename =  strcat('Params/simulationGibbs_',data,'_Simple');
%saveas(gcf,figurename,'epsc')
