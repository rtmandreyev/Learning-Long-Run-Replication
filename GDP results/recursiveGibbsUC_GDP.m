function [muDraws,gamDraws,rho1Draws,rho2Draws,sigDraws,zRealTimeDraws,...
    xRealTimeDraws,zSmoothedDraws,xSmoothedDraws,yForecasts,yGrowthForecastsA] = ...
    recursiveGibbsUC_GDP(y,muParams,gamParams,rhoParams,rho2Params,sig2Params,...
    burnIn,B,estimationStart,vintageStart,stepSize,paramInd,seed,fullFlag)

rng(seed,'Philox')

T = size(y,1);
muDraws = NaN(T,B);
rho1Draws = NaN(T,B);
rho2Draws = NaN(T,B);
gamDraws = NaN(T,B);
sigDraws = NaN(T,B);
zRealTimeDraws = NaN(T,B);
xRealTimeDraws = NaN(T,B);
xLRealTimeDraws = NaN(T,B);
sigGamPropInit = 0.05;

maxHorizon = 4*11+1;
zForecasts = NaN(T,maxHorizon,B);
xForecasts = NaN(T,maxHorizon,B);
xLForecasts = NaN(T,maxHorizon,B);
curVintage = vintageStart;
startInd = 51;

initVar = 0.01^2;
initVarMat = [initVar,-initVar,-0.9*initVar;-initVar,initVar,0.9*initVar;...
    -0.9*initVar,0.9*initVar,initVar];

for t = estimationStart:T
    
    t
    
    if size(y,2) > 1
        yt = y(startInd:t,curVintage);
    else
        yt = y(startInd:t);
    end
    
    if ~mod(t-estimationStart,stepSize) && paramInd(t)
        
        if t > estimationStart
            [muNext,gamNext,rho1Next,rho2Next,sigNext,zNext,xNext,xLNext,sigGamProp] = ...
                gibbsUC_GDP(yt,muParams,gamParams,rhoParams,rho2Params,...
                sig2Params,burnIn,B,sigGamPropInit,paramInd(startInd:t),...
                meanNext);
        else
            [muNext,gamNext,rho1Next,rho2Next,sigNext,zNext,xNext,xLNext,sigGamProp] = ...
                gibbsUC_GDP(yt,muParams,gamParams,rhoParams,rho2Params,...
                sig2Params,burnIn,B,sigGamPropInit,paramInd(startInd:t),[]);
        end
        meanNext = [mean(muNext),mean(gamNext),mean(rho1Next),...
                mean(rho2Next),mean(sigNext.^2)];
        useInd = (abs((rho1Next+sqrt(rho1Next.^2+4*rho2Next))./2) < 1) & ...
            (abs((rho1Next-sqrt(rho1Next.^2+4*rho2Next))./2) < 1) & ...
            (muNext > 0);
        sum(useInd)
        
        figure(1)
        tiledlayout(3,2)
        nexttile
        histogram(muNext)
        nexttile
        histogram(rho1Next);
        nexttile
        histogram(rho2Next)
        nexttile
        histogram(gamNext)
        nexttile
        histogram(sigNext)
        
        muDraws(t,useInd) = muNext(useInd);
        gamDraws(t,useInd) = gamNext(useInd);
        rho1Draws(t,useInd) = rho1Next(useInd);
        rho2Draws(t,useInd) = rho2Next(useInd);
        sigDraws(t,useInd) = sigNext(useInd);
        zRealTimeDraws(t,useInd) = zNext(useInd,end);
        xRealTimeDraws(t,useInd) = xNext(useInd,end);
        xLRealTimeDraws(t,useInd) = xLNext(useInd,end);

        sigGamPropInit = sigGamProp;
        
    else
        
        muDraws(t,:) = muDraws(t-1,:);
        gamDraws(t,:) = gamDraws(t-1,:);
        rho1Draws(t,:) = rho1Draws(t-1,:);
        rho2Draws(t,:) = rho2Draws(t-1,:);
        sigDraws(t,:) = sigDraws(t-1,:);
        zNext = NaN(B,t);
        xNext = NaN(B,t);
        xLNext = NaN(B,t);
        if fullFlag || (t==T)
            for b = 1:B
                if ~isnan(rho1Draws(t,b))
                    F = [1 0 0 1; 0 rho1Draws(t,b) rho2Draws(t,b) 0; 0 1 0 0; 0 0 0 1];
                    Q = [sqrt(gamDraws(t,b))*sigDraws(t,b) 0 0 0;...
                        0 sqrt((1-gamDraws(t,b)))*sigDraws(t,b) 0 0;...
                        0 0 0 0; 0 0 0 0];
                    H = [1 1 0 0];
                    Mdl = ssm(F,Q,H,'Mean0',[yt(1);0;0;muDraws(t,b)],'Cov0',...
                        [[initVarMat,zeros(3,1)];zeros(1,4)]);
                    smoothedDraws = simsmooth(Mdl,yt);
                    zNext(b,startInd:t) = smoothedDraws(:,1)';
                    xNext(b,startInd:t) = smoothedDraws(:,2)';
                    xLNext(b,startInd:t) = smoothedDraws(:,3)';
                end
            end
            zRealTimeDraws(t,:) = zNext(:,end);
            xRealTimeDraws(t,:) = xNext(:,end);
            xLRealTimeDraws(t,:) = xLNext(:,end);
        end
        
    end
    
    zForecasts(t,1,:) = zNext(:,end);
    xForecasts(t,1,:) = xNext(:,end);
    xLForecasts(t,1,:) = xLNext(:,end);
    forecastShocks1 = sqrt(gamDraws(t,:)).*sigDraws(t,:).*randn(maxHorizon,B);
    forecastShocks2 = sqrt(1-gamDraws(t,:)).*sigDraws(t,:).*randn(maxHorizon,B);
    for ii = 1:(maxHorizon-1)
        
        zForecasts(t,ii+1,:) = muDraws(t,:)' + squeeze(zForecasts(t,ii,:)) + forecastShocks1(ii+1,:)';
        xForecasts(t,ii+1,:) = rho1Draws(t,:)'.*squeeze(xForecasts(t,ii,:)) + ...
            rho2Draws(t,:)'.*squeeze(xLForecasts(t,ii,:)) + forecastShocks2(ii+1,:)';
        xLForecasts(t,ii+1,:) = xForecasts(t,ii,:);
        
    end
    
    curVintage = curVintage + 1;
    
    if t > T-stepSize
        
        zSmoothedDraws = [NaN(startInd-1,B); zNext'];
        xSmoothedDraws = [NaN(startInd-1,B); xNext'];
        
    end
    
end

yForecasts = zForecasts + xForecasts;

yGrowthForecastsA = NaN(size(yForecasts(estimationStart:4:end,:,:),1),11,B);
curInd = 1;
curVintage = vintageStart;
for t = estimationStart:4:T
    
    yt = y(1:t,curVintage);
    yLevelCur = mean(exp(yt(t-3:t)));
    for ii = 1:11
        
        yLevelNext = squeeze(mean(exp(yForecasts(t,4*(ii-1)+2:4*ii+1,:)),2));
        yGrowthForecastsA(curInd,ii,:) = 100*(yLevelNext./yLevelCur-1);
        yLevelCur = yLevelNext;
        
    end
    curVintage = curVintage + 4;
    curInd = curInd + 1;
end

end