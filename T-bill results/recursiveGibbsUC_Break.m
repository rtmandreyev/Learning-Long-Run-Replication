function [rhoDraws,gamDraws,sigDraws,muRealTimeDraws,xRealTimeDraws,...
    muSmoothedDraws,xSmoothedDraws,yForecasts,modelYields] = ...
    recursiveGibbsUC_Break(y,rhoParams,gam1Params,gam2Params,sig2Params,...
    burnIn,B,estimationStart,stepSize,paramInd,seed,breakInd)

rng(seed,'Philox')

T = numel(y);
rhoDraws = NaN(T,B);
gamDraws = NaN(T,B);
sigDraws = NaN(T,B);
muRealTimeDraws = NaN(T,B);
xRealTimeDraws = NaN(T,B);
sigGamPropInit = 0.1;

maxHorizon = 4*30;
muForecasts = NaN(T,maxHorizon,B);
xForecasts = NaN(T,maxHorizon,B);
modelYields = NaN(T,120);

for t = estimationStart:T
    
    t
    
    if ~mod(t-estimationStart,stepSize) && paramInd(t)
        
        if (t > estimationStart) && (t < breakInd)
            [rhoNext,gamNext,sigNext,muNext,xNext,sigGamProp] = gibbsUCSmooth(y(1:t),rhoParams,...
                gam1Params,sig2Params,burnIn,B,sigGamPropInit,paramInd(1:t),...
                [mean(rhoDraws(t-1,:),'omitnan'),mean(gamDraws(t-1,:),'omitnan'),...
                mean(sigDraws(t-1,:).^2,'omitnan')]);
        elseif (t >= breakInd)
            [rhoNext,gamNext,sigNext,muNext,xNext,sigGamProp] = gibbsUCSmooth(y(1:t),rhoParams,...
                gam2Params,sig2Params,burnIn,B,sigGamPropInit,paramInd(1:t),...
                [mean(rhoDraws(t-1,:),'omitnan'),mean(gamDraws(t-1,:),'omitnan'),...
                mean(sigDraws(t-1,:).^2,'omitnan')]);
        else
            [rhoNext,gamNext,sigNext,muNext,xNext,sigGamProp] = gibbsUCSmooth(y(1:t),rhoParams,...
                gam1Params,sig2Params,burnIn,B,sigGamPropInit,paramInd(1:t),[]);
        end
        stationaryInd = (rhoNext<1);
        rhoDraws(t,stationaryInd) = rhoNext(stationaryInd);
        gamDraws(t,stationaryInd) = gamNext(stationaryInd);
        sigDraws(t,stationaryInd) = sigNext(stationaryInd);
        muDrawsCur = muNext;
        xDrawsCur = xNext;
        sigGamPropInit = sigGamProp;

        figure(1)
        tiledlayout(3,1)
        nexttile
        histogram(rhoNext)
        axis tight
        nexttile
        histogram(gamNext)
        axis tight
        nexttile
        histogram(sigNext)
        axis tight
        
    else
        
        rhoDraws(t,:) = rhoDraws(t-1,:);
        gamDraws(t,:) = gamDraws(t-1,:);
        sigDraws(t,:) = sigDraws(t-1,:);
        
        muDrawsNext = NaN(size(muDrawsCur,1),size(muDrawsCur,2)+1);
        xDrawsNext = NaN(size(xDrawsCur,1),size(xDrawsCur,2)+1);
        muDrawsNext(:,1) = muDrawsCur(:,1);
        xDrawsNext(:,1) = xDrawsCur(:,1);
        muDrawsCur = muDrawsNext;
        xDrawsCur = xDrawsNext;
        
        for b = 1:B
            if ~isnan(rhoDraws(t,b))
                F = [1 0; 0 rhoDraws(t,b)];
                Q = [sqrt(gamDraws(t,b)*sigDraws(t,b)^2) 0; 0 sqrt((1-gamDraws(t,b))*sigDraws(t,b)^2)];
                H = [1 1];
                Mdl = ssm(F,Q,H,'Mean0',[y(1);0],'Cov0',[1,-1;-1,1]);
                smoothedDraws = simsmooth(Mdl,y(1:t));
                muDrawsCur(b,:) = smoothedDraws(:,1);
                xDrawsCur(b,:) = smoothedDraws(:,2);
            end
        end
        
    end
    
    muRealTimeDraws(t,:) = muDrawsCur(:,t);
    xRealTimeDraws(t,:) = xDrawsCur(:,t);
    
    if t == T
        
        muSmoothedDraws = muDrawsCur';
        xSmoothedDraws = xDrawsCur';
        
    end
    
    muForecasts(t,1,:) = muDrawsCur(:,t);
    xForecasts(t,1,:) = xDrawsCur(:,t);
    forecastShocks1 = sqrt(gamDraws(t,:)).*sigDraws(t,:).*randn(maxHorizon,B);
    forecastShocks2 = sqrt(1-gamDraws(t,:)).*sigDraws(t,:).*randn(maxHorizon,B);

    for ii = 1:(maxHorizon-1)
        
        muForecasts(t,ii+1,:) = squeeze(muForecasts(t,ii,:)) + forecastShocks1(ii+1,:)';
        xForecasts(t,ii+1,:) = rhoDraws(t,:)'.*squeeze(xForecasts(t,ii,:)) + forecastShocks2(ii+1,:)';
        
    end
    
end

for t = estimationStart:T
    t
    muTemp = [squeeze(muForecasts(t,:,:))',NaN(B,80)];
    xTemp = [squeeze(xForecasts(t,:,:))',NaN(B,80)];
    for ii = maxHorizon:119
        muTemp(:,ii+1) = muTemp(:,ii) + sqrt(gamDraws(t,:)').*randn(B,1);
        xTemp(:,ii+1) = rhoDraws(t,:)'.*xTemp(:,ii) + sqrt(1-gamDraws(t,:)').*randn(B,1);
        %muForecasts(t,ii+1,:) = squeeze(muForecasts(t,ii,:)) + sqrt(gamDraws(t,:)').*randn(B,1);
        %xForecasts(t,ii+1,:) = rhoDraws(t,:)'.*squeeze(xForecasts(t,ii,:)) + sqrt(1-gamDraws(t,:)').*randn(B,1);
    end
    modelYields(t,:) = cumsum(mean(muTemp+xTemp,'omitnan'),2)./(1:120);
end

yForecasts = muForecasts + xForecasts;
%modelYields = cumsum(squeeze(mean(yForecasts,3,'omitnan')),2)./(1:120);
%yForecasts = yForecasts(:,1:40,:);

end