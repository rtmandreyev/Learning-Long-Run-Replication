function [muDraws,gamDraws,rho1Draws,rho2Draws,sigDraws,zDraws,xDraws,xLDraws,sigGamProp] = ...
    gibbsUC_GDP(y,muParams,gamParams,rhoParams,rho2Params,...
    sig2Params,burnIn,B,sigGamPropInit,paramInd,initDraws)

T = numel(y);
T1 = sum(paramInd);
paramInd = logical(paramInd);

muMu = muParams(1);
varMu = muParams(2);

muRho = [rhoParams(1)-rho2Params(1); rho2Params(1)];
varRho = [rhoParams(2)+rho2Params(2),-rho2Params(2); -rho2Params(2),rho2Params(2)];

aGam = gamParams(1);
bGam = gamParams(2);
logBetaPdf = @(x,a,b) (a-1)*log(x) + (b-1)*log(1-x) - betaln(a,b);
logGamPrior = @(x) logBetaPdf(x,aGam,bGam);
gamAcceptance = zeros(burnIn,1);
sigGamProp = sigGamPropInit;

aSig2 = sig2Params(1);
bSig2 = sig2Params(2);

if isempty(initDraws)
    muDraws = [muMu; NaN(burnIn+B-1,1)];
    rho1Draws = [muRho(1); NaN(burnIn+B-1,1)];
    rho2Draws = [muRho(2); NaN(burnIn+B-1,1)];
    gamDraws = [(aGam-1)/(aGam+bGam-2) ;NaN(burnIn+B-1,1)];
    sig2Draws = [bSig2/(aSig2+1); NaN(burnIn+B-1,1)];
else
    muDraws = [initDraws(1); NaN(burnIn+B-1,1)];
    rho1Draws = [initDraws(3); NaN(burnIn+B-1,1)];
    rho2Draws = [initDraws(4); NaN(burnIn+B-1,1)];
    gamDraws = [initDraws(2);NaN(burnIn+B-1,1)];
    sig2Draws = [initDraws(5); NaN(burnIn+B-1,1)];
end

F = [1 0 0 1; 0 rho1Draws(1) rho2Draws(1) 0; 0 1 0 0; 0 0 0 1];
Q = [sqrt(gamDraws(1)*sig2Draws(1)) 0 0 0; 0 sqrt((1-gamDraws(1))*sig2Draws(1)) 0 0;...
    0 0 0 0; 0 0 0 0];
H = [1 1 0 0];
initVar = 0.01^2;
initVarMat = [initVar,-initVar,-0.9*initVar;-initVar,initVar,0.9*initVar;...
    -0.9*initVar,0.9*initVar,initVar];
Mdl = ssm(F,Q,H,'Mean0',[y(1);0;0;muDraws(1)],'Cov0',[[initVarMat,zeros(3,1)];zeros(1,4)]);

statesSmoothedGuess = simsmooth(Mdl,y);
zDraws = [statesSmoothedGuess(:,1)'; NaN(burnIn+B-1,T)];
xDraws = [statesSmoothedGuess(:,2)'; NaN(burnIn+B-1,T)];
xLDraws = [statesSmoothedGuess(:,3)'; NaN(burnIn+B-1,T)];

paramInd1 = logical([paramInd(1:end-2);0;0]);
paramInd2 = logical([0;paramInd(2:end-1);0]);
paramInd3 = logical([0;0;paramInd(3:end)]);

while true
    
    for b = (1:burnIn-1)
        
        % Draw mu
        
        varMuPosterior = 1/(1/varMu + (T1-1)/(gamDraws(b)*sig2Draws(b)));
        muMuPosterior = varMuPosterior*(muMu/varMu + sum(diff(zDraws(b,paramInd)))/...
            (gamDraws(b)*sig2Draws(b)));
        muDraws(b+1) = normrnd(muMuPosterior,sqrt(varMuPosterior));
        
        % Draw rho
        
        varRhoPosterior = pinv(pinv(varRho) + [xDraws(b,paramInd2);...
            xDraws(b,paramInd1)]*[xDraws(b,paramInd2);xDraws(b,paramInd1)]'./...
            ((1-gamDraws(b))*sig2Draws(b)));
        muRhoPosterior = varRhoPosterior*(pinv(varRho)*muRho + [xDraws(b,paramInd2);...
            xDraws(b,paramInd1)]*xDraws(b,paramInd3)'./((1-gamDraws(b))*sig2Draws(b)));
        corrRho = varRhoPosterior(1,2)/sqrt(varRhoPosterior(1,1)*varRhoPosterior(2,2));
        rhoDrawsCur(1) = normrnd(muRhoPosterior(1),sqrt(varRhoPosterior(1,1)));
        rhoDrawsCur(2) = normrnd(muRhoPosterior(2) + ...
            sqrt(varRhoPosterior(2,2))/sqrt(varRhoPosterior(1,1))*...
            corrRho*(rhoDrawsCur(1)-muRhoPosterior(1)),sqrt((1-corrRho^2)*varRhoPosterior(2,2)));
        if (rhoDrawsCur(1)+rhoDrawsCur(2)<1) && (abs((rhoDrawsCur(1)+...
                sqrt(rhoDrawsCur(1).^2+4*rhoDrawsCur(2)))./2) < 1) && ...
                (abs((rhoDrawsCur(1)-sqrt(rhoDrawsCur(1).^2+4*rhoDrawsCur(2)))./2) < 1)
            rho1Draws(b+1) = rhoDrawsCur(1);
            rho2Draws(b+1) = rhoDrawsCur(2);
        else
            rho1Draws(b+1) = rho1Draws(b);
            rho2Draws(b+1) = rho2Draws(b);
        end
        
        % Draw gamma
        
        gamProp = normrnd(gamDraws(b),sigGamProp);
        
        logLCur = -(T1-1)/2*log(sig2Draws(b)) - (T1-1)/2*log(gamDraws(b)) - ...
            sum(((diff(zDraws(b,paramInd))-muDraws(b+1))./sqrt(sig2Draws(b)*gamDraws(b))).^2)/2 - ...
            (T1-2)/2*log(sig2Draws(b)) - (T1-2)/2*log(1-gamDraws(b)) - ...
            sum(((xDraws(b,paramInd3)-rho1Draws(b+1)*xDraws(b,paramInd2)-...
            rho2Draws(b+1)*xDraws(b,paramInd1))./sqrt(sig2Draws(b)*(1-gamDraws(b)))).^2)/2 + ...
            logGamPrior(gamDraws(b));
        logLProp = -(T1-1)/2*log(sig2Draws(b)) - (T1-1)/2*log(gamProp) - ...
            sum(((diff(zDraws(b,paramInd))-muDraws(b+1))./sqrt(sig2Draws(b)*gamProp)).^2)/2 - ...
            (T1-2)/2*log(sig2Draws(b)) - (T1-2)/2*log(1-gamProp) - ...
            sum(((xDraws(b,paramInd3)-rho1Draws(b+1)*xDraws(b,paramInd2)-...
            rho2Draws(b+1)*xDraws(b,paramInd1))./sqrt(sig2Draws(b)*(1-gamProp))).^2)/2 + ...
            logGamPrior(gamProp);
        
        if isreal(logLProp)
            alph = exp(logLProp-logLCur);
        else
            alph = 0;
        end
        if rand <= alph
            gamDraws(b+1) = gamProp;
            gamAcceptance(b) = 1;
        else
            gamDraws(b+1) = gamDraws(b);
        end
        
        % Draw sigma^2
        
        aSig2Posterior = aSig2 + (2*T1-3)/2;
        bSig2Posterior = bSig2 + sum(((diff(zDraws(b,paramInd))-muDraws(b+1))./sqrt(gamDraws(b+1))).^2)/2 + ...
            sum(((xDraws(b,paramInd3)-rho1Draws(b+1)*xDraws(b,paramInd2)-...
            rho2Draws(b+1)*xDraws(b,paramInd1))./sqrt(1-gamDraws(b+1))).^2)/2;
        sig2Draws(b+1) = 1/gamrnd(aSig2Posterior,1/bSig2Posterior);
        
        % Draw mu and x
        
        F = [1 0 0 1; 0 rho1Draws(b+1) rho2Draws(b+1) 0; 0 1 0 0; 0 0 0 1];
        Q = [sqrt(gamDraws(b+1)*sig2Draws(b+1)) 0 0 0; 0 sqrt((1-gamDraws(b+1))*sig2Draws(b+1)) 0 0;...
            0 0 0 0; 0 0 0 0];
        H = [1 1 0 0];
        Mdl = ssm(F,Q,H,'Mean0',[y(1);0;0;muDraws(b+1)],'Cov0',[[initVarMat,zeros(3,1)];zeros(1,4)]);
        smoothedDraws = simsmooth(Mdl,y);
        zDraws(b+1,:) = smoothedDraws(:,1)';
        xDraws(b+1,:) = smoothedDraws(:,2)';
        xLDraws(b+1,:) = smoothedDraws(:,3)';
        
    end
    
    mean(gamAcceptance)
    if mean(gamAcceptance) < 0.25
        sigGamProp = 0.8*sigGamProp;
        gamAcceptance = zeros(burnIn,1);
    elseif mean(gamAcceptance) > 0.4
        sigGamProp = 1.2*sigGamProp;
        gamAcceptance = zeros(burnIn,1);
    else
        break
    end
    muDraws(1) = muDraws(burnIn);
    rho1Draws(1) = rho1Draws(burnIn);
    rho2Draws(1) = rho2Draws(burnIn);
    gamDraws(1) = gamDraws(burnIn);
    sig2Draws(1) = sig2Draws(burnIn);
    zDraws(1,:) = zDraws(burnIn,:);
    xDraws(1,:) = xDraws(burnIn,:);
    xLDraws(1,:) = xLDraws(burnIn,:);
    
end

for b = (burnIn:burnIn+B-1)
    
    % Draw mu
    
    varMuPosterior = 1/(1/varMu + (T1-1)/(gamDraws(b)*sig2Draws(b)));
    muMuPosterior = varMuPosterior*(muMu/varMu + sum(diff(zDraws(b,paramInd)))/...
        (gamDraws(b)*sig2Draws(b)));
    muDraws(b+1) = normrnd(muMuPosterior,sqrt(varMuPosterior));
    
    % Draw rho
    
    varRhoPosterior = pinv(pinv(varRho) + [xDraws(b,paramInd2);...
        xDraws(b,paramInd1)]*[xDraws(b,paramInd2);xDraws(b,paramInd1)]'./...
        ((1-gamDraws(b))*sig2Draws(b)));
    muRhoPosterior = varRhoPosterior*(pinv(varRho)*muRho + [xDraws(b,paramInd2);...
        xDraws(b,paramInd1)]*xDraws(b,paramInd3)'./((1-gamDraws(b))*sig2Draws(b)));
    corrRho = varRhoPosterior(1,2)/sqrt(varRhoPosterior(1,1)*varRhoPosterior(2,2));
    rhoDrawsCur(1) = normrnd(muRhoPosterior(1),sqrt(varRhoPosterior(1,1)));
    rhoDrawsCur(2) = normrnd(muRhoPosterior(2) + ...
        sqrt(varRhoPosterior(2,2))/sqrt(varRhoPosterior(1,1))*...
        corrRho*(rhoDrawsCur(1)-muRhoPosterior(1)),sqrt((1-corrRho^2)*varRhoPosterior(2,2)));
    if (rhoDrawsCur(1)+rhoDrawsCur(2)<1) && (abs((rhoDrawsCur(1)+...
                sqrt(rhoDrawsCur(1).^2+4*rhoDrawsCur(2)))./2) < 1) && ...
                (abs((rhoDrawsCur(1)-sqrt(rhoDrawsCur(1).^2+4*rhoDrawsCur(2)))./2) < 1)
        rho1Draws(b+1) = rhoDrawsCur(1);
        rho2Draws(b+1) = rhoDrawsCur(2);
    else
        rho1Draws(b+1) = rho1Draws(b);
        rho2Draws(b+1) = rho2Draws(b);
    end
    
    % Draw gamma
    
    gamProp = normrnd(gamDraws(b),sigGamProp);
    
    logLCur = -(T1-1)/2*log(sig2Draws(b)) - (T1-1)/2*log(gamDraws(b)) - ...
        sum(((diff(zDraws(b,paramInd))-muDraws(b+1))./sqrt(sig2Draws(b)*gamDraws(b))).^2)/2 - ...
        (T1-2)/2*log(sig2Draws(b)) - (T1-2)/2*log(1-gamDraws(b)) - ...
        sum(((xDraws(b,paramInd3)-rho1Draws(b+1)*xDraws(b,paramInd2)-...
        rho2Draws(b+1)*xDraws(b,paramInd1))./sqrt(sig2Draws(b)*(1-gamDraws(b)))).^2)/2 + ...
        logGamPrior(gamDraws(b));
    logLProp = -(T1-1)/2*log(sig2Draws(b)) - (T1-1)/2*log(gamProp) - ...
        sum(((diff(zDraws(b,paramInd))-muDraws(b+1))./sqrt(sig2Draws(b)*gamProp)).^2)/2 - ...
        (T1-2)/2*log(sig2Draws(b)) - (T1-2)/2*log(1-gamProp) - ...
        sum(((xDraws(b,paramInd3)-rho1Draws(b+1)*xDraws(b,paramInd2)-...
        rho2Draws(b+1)*xDraws(b,paramInd1))./sqrt(sig2Draws(b)*(1-gamProp))).^2)/2 + ...
        logGamPrior(gamProp);
    
    if isreal(logLProp)
        alph = exp(logLProp-logLCur);
    else
        alph = 0;
    end
    if rand <= alph
        gamDraws(b+1) = gamProp;
        gamAcceptance(b) = 1;
    else
        gamDraws(b+1) = gamDraws(b);
    end
    
    % Draw sigma^2
    
    aSig2Posterior = aSig2 + (2*T1-3)/2;
    bSig2Posterior = bSig2 + sum(((diff(zDraws(b,paramInd))-muDraws(b+1))./sqrt(gamDraws(b+1))).^2)/2 + ...
        sum(((xDraws(b,paramInd3)-rho1Draws(b+1)*xDraws(b,paramInd2)-...
        rho2Draws(b+1)*xDraws(b,paramInd1))./sqrt(1-gamDraws(b+1))).^2)/2;
    sig2Draws(b+1) = 1/gamrnd(aSig2Posterior,1/bSig2Posterior);
    
    % Draw mu and x
    
    F = [1 0 0 1; 0 rho1Draws(b+1) rho2Draws(b+1) 0; 0 1 0 0; 0 0 0 1];
    Q = [sqrt(gamDraws(b+1)*sig2Draws(b+1)) 0 0 0; 0 sqrt((1-gamDraws(b+1))*sig2Draws(b+1)) 0 0;...
        0 0 0 0; 0 0 0 0];
    H = [1 1 0 0];
    Mdl = ssm(F,Q,H,'Mean0',[y(1);0;0;muDraws(b+1)],'Cov0',[[initVarMat,zeros(3,1)];zeros(1,4)]);
    smoothedDraws = simsmooth(Mdl,y);
    zDraws(b+1,:) = smoothedDraws(:,1)';
    xDraws(b+1,:) = smoothedDraws(:,2)';
    xLDraws(b+1,:) = smoothedDraws(:,3)';
    
end

muDraws = muDraws(burnIn+1:end,:);
gamDraws = gamDraws(burnIn+1:end);
rho1Draws = rho1Draws(burnIn+1:end);
rho2Draws = rho2Draws(burnIn+1:end);
sig2Draws = sig2Draws(burnIn+1:end);
zDraws = zDraws(burnIn+1:end,:);
xDraws = xDraws(burnIn+1:end,:);
xLDraws = xLDraws(burnIn+1:end,:);

sigDraws = sqrt(sig2Draws);