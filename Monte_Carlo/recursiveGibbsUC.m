function [rhoDraws,gamDraws,sigDraws,muRealTimeDraws,xRealTimeDraws,...
    muSmoothedDraws,xSmoothedDraws,yForecasts] = ...
    recursiveGibbsUC(y,rhoParams,gamParams,sig2Params,burnIn,numKeep,estimationStart,...
    stepSize,paramInd,seed)

rng(seed,'Philox')

T = numel(y);
rhoDraws = NaN(T,numKeep);
gamDraws = NaN(T,numKeep);
sigDraws = NaN(T,numKeep);
muRealTimeDraws = NaN(T,numKeep);
xRealTimeDraws = NaN(T,numKeep);
sigGamPropInit = 0.1;

maxHorizon = 4*10;
muForecasts = NaN(T,maxHorizon,numKeep);
xForecasts = NaN(T,maxHorizon,numKeep);

for t = estimationStart:T
    
    %t
    
    if ~mod(t-estimationStart,stepSize) && paramInd(t)
        
        if t > estimationStart
            [rhoNext,gamNext,sigNext,muNext,xNext,sigGamProp] = gibbsUC(y(1:t),rhoParams,...
                gamParams,sig2Params,burnIn,numKeep,sigGamPropInit,paramInd(1:t),...
                [mean(rhoDraws(t-1,:),'omitnan'),mean(gamDraws(t-1,:),'omitnan'),...
                mean(sigDraws(t-1,:),'omitnan')]);
        else
            [rhoNext,gamNext,sigNext,muNext,xNext,sigGamProp] = gibbsUC(y(1:t),rhoParams,...
                gamParams,sig2Params,burnIn,numKeep,sigGamPropInit,paramInd(1:t),[]);
        end
        
        stationaryInd = (rhoNext<1);
        rhoDraws(t,stationaryInd) = rhoNext(stationaryInd);
        gamDraws(t,stationaryInd) = gamNext(stationaryInd);
        sigDraws(t,stationaryInd) = sigNext(stationaryInd);
        muDrawsCur = muNext;
        xDrawsCur = xNext;
        sigGamPropInit = sigGamProp;
        
%         close all;
%         figure(1)
%         tiledlayout(3,1)
%         nexttile
%         histogram(rhoNext)
%         axis tight
%         nexttile
%         histogram(gamNext)
%         axis tight
%         nexttile
%         histogram(sigNext)
%         axis tight
        
        rhoNext;
        
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
        
        for b = 1:numKeep
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
    forecastShocks1 = sqrt(gamDraws(t,:)).*sigDraws(t,:).*randn(maxHorizon,numKeep);
    forecastShocks2 = sqrt(1-gamDraws(t,:)).*sigDraws(t,:).*randn(maxHorizon,numKeep);
    for ii = 1:(maxHorizon-1)
        
        muForecasts(t,ii+1,:) = squeeze(muForecasts(t,ii,:)) + forecastShocks1(ii+1,:)';
        xForecasts(t,ii+1,:) = rhoDraws(t,:)'.*squeeze(xForecasts(t,ii,:)) + forecastShocks2(ii+1,:)';
        
    end
    
end

yForecasts = muForecasts + xForecasts;

end
