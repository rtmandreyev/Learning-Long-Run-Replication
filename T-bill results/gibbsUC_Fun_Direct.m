function [obj,yForecasts,rhoDraws,gamDraws,meanForecasts] = gibbsUC_Fun_Direct(thet,seed)

try
    
    load('Tbl_Data.mat','tblData','spfForecastsQ','T','yieldData','spfForecastsY',...
        'spfForecasts10','dateVec','estimationStart')
    load('Tbl_Regression_Stats.mat','biasMat','biasMat_se','arMat','arMat_se',...
        'mzMat','mzMat_se','cgMat','cgMat_se','csConventionalMat',...
        'csConventionalMat_se','csContrarianMat','csContrarianMat_se')
    
    paramInd = (tblData>0.25);
    rhoParams = thet(1:2);
    gamMoments = thet(3:4);
    gamFun = @(x) [(x(1)-1)/(sum(x)-2); sqrt(prod(x)/((sum(x)^2)*(sum(x)+1)))];
    options = optimoptions(@fminunc,'Display','off');
    if size(gamMoments,1) < size(gamMoments,2)
        gamMoments = gamMoments';
    end
    momentFun = @(x) 1000*sum(((gamFun(x)-gamMoments)./gamMoments).^2);
    gamParams = fminunc(momentFun,[2;2],options);
    if any(isnan(gamParams))
        obj = NaN
        return
    end
    sig2Params = [1.25;0.5625];
    
    burnIn = 25000;
    B = 25000;
    stepSize = 4;
    yieldInd = [2 3 4 8 12 20 40];
    
    [~,~,~,~,~,~,~,yForecasts] = recursiveGibbsUC(tblData,rhoParams,gamParams,...
        sig2Params,burnIn,B,estimationStart,stepSize,paramInd,seed);
    meanForecasts = squeeze(mean(yForecasts,3,'omitnan'));
    modelYields = cumsum(meanForecasts,2)./(1:40);
    
    if sum(all(isnan(meanForecasts),2)) == T
        obj = NaN
        return
    end
    
    %     modelForecastsQ = meanForecasts(:,2:5);
    %
    %     modelForecastsY = NaN(size(spfForecastsY));
    %     modelForecasts10 = NaN(size(spfForecasts10));
    %     for t = estimationStart:T
    %         switch quarter(dateVec(t))
    %             case 1
    %                 modelForecastsY(t,:) = [mean([tblData(t),meanForecasts(t,1:3)]),...
    %                     mean(meanForecasts(t,4:7)),mean(meanForecasts(t,8:11)),...
    %                     mean(meanForecasts(t,12:15))];
    %                 modelForecasts10(t) = mean([tblData(t),meanForecasts(t,1:39)]);
    %             case 2
    %                 modelForecastsY(t,:) = [mean([tblData(t-1:t)',meanForecasts(t,1:2)]),...
    %                     mean(meanForecasts(t,3:6)),mean(meanForecasts(t,7:10)),...
    %                     mean(meanForecasts(t,11:14))];
    %             case 3
    %                 modelForecastsY(t,:) = [mean([tblData(t-2:t)',meanForecasts(t,1)]),...
    %                     mean(meanForecasts(t,2:5)),mean(meanForecasts(t,6:9)),...
    %                     mean(meanForecasts(t,10:13))];
    %             case 4
    %                 modelForecastsY(t,:) = [mean(tblData(t-3:t)),...
    %                     mean(meanForecasts(t,1:4)),mean(meanForecasts(t,5:8)),...
    %                     mean(meanForecasts(t,9:12))];
    %         end
    %     end
    
    %     modelYields(isnan(yieldData(:,1:40))) = NaN;
    %     obj = sum(sum((modelForecastsQ-spfForecastsQ(:,2:end)).^2,'omitnan')) + ...
    %         sum(sum((modelForecastsY-spfForecastsY).^2,'omitnan')) + ...
    %         sum((modelForecasts10-spfForecasts10).^2,'omitnan') + ...
    %         sum(sum(((yieldData(:,yieldInd)-mean(yieldData(:,yieldInd),'omitnan'))...
    %         -(modelYields(:,yieldInd)-mean(modelYields(:,yieldInd),'omitnan'))).^2,'omitnan'))
    
    maxHorizon = 7;
    spfInd = ~isnan(spfForecastsQ(:,1));
    meanForecasts(~spfInd,:) = NaN;
    staggeredModelForecasts = meanForecasts(:,1:(maxHorizon+1));
    
    for ii = 0:maxHorizon
        
        staggeredModelForecasts(:,ii+1) = circshift(meanForecasts(:,ii+1),ii);
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
    
    %[biasMatModel, biasMatModel_se, biasMatModel_p]'
    %[arMatModel(:,2), arMatModel_se(:,2), arMatModel_p(:,2)]'
    %[mzMatModel(:,2), mzMatModel_se(:,2), mzMatModel_p(:,2)]'
    %[cgMatModel(:,2), cgMatModel_se(:,2), cgMatModel_p(:,2)]'
    
    %% Compute expectations hypothesis violations
    
    maxHorizon = 40;
    
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
    
    %[csConventionalMatModel; csConventionalMatModel_se; csConventionalMatModel_p]
    %[csContrarianMatModel; csContrarianMatModel_se; csContrarianMatModel_p]
    tempAR = squeeze(arMat(2:end,:,1)) - arMatModel(2:end,:);
    tempARRef = squeeze(arMat(2:end,:,1));
    tempAR_se = squeeze(arMat_se(2:end,:,1)) - arMatModel_se(2:end,:);
    tempMZ = squeeze(mzMat(2:end,:,1)) - mzMatModel(2:end,:);
    tempMZRef = squeeze(mzMat(2:end,:,1));
    tempMZ_se = squeeze(mzMat_se(2:end,:,1)) - mzMatModel_se(2:end,:);
    tempCG = squeeze(cgMat(2:end,:,1)) - cgMatModel(2:end,:);
    tempCGRef = squeeze(cgMat(2:end,:,1));
    tempCG_se = squeeze(cgMat_se(2:end,:,1)) - cgMatModel_se(2:end,:);
    
    %     obj = sum([(biasMat(2:end,1)-biasMatModel(2:end))./biasMat(2:end,1);...
    %         tempAR(:)./tempARRef(:);...
    %         tempMZ(:)./tempMZRef(:);...
    %         tempCG(:)./tempCGRef(:);...
    %         (csConventionalMat'-csConventionalMatModel')./(csConventionalMat');...
    %         (csContrarianMat'-csContrarianMatModel')./(csContrarianMat')].^2,...
    %         'omitnan')

    obj = sum(sum([meanForecasts(:,2:5)-spfForecastsQ(:,2:5),...
        modelYields(:,[20,40])-yieldData(:,[20,40])].^2,'omitnan'))
    
    if obj == 0
        obj = NaN
    end
    
catch
    
    obj = NaN
    
end

end