function [rhoDraws,gamDraws,sigDraws,muDraws,xDraws,sigGamProp] = ...
    gibbsUCSmooth(y,rhoParams,gamParams,sig2Params,burnIn,B,...
    sigGamPropInit,paramInd,initDraws)

T = numel(y);
numInd = find(paramInd);
keepInd = (numInd(2:end)-numInd(1:end-1)) == 1;
T1 = sum(paramInd) - sum(keepInd==0);
paramInd = logical(paramInd);

muRho = rhoParams(1);
varRho = rhoParams(2);

aGam = gamParams(1);
bGam = gamParams(2);
logBetaPdf = @(x,a,b) (a-1)*log(x) + (b-1)*log(1-x) - betaln(a,b);
logGamPrior = @(x) logBetaPdf(x,aGam,bGam);
gamAcceptance = zeros(burnIn,1);
sigGamProp = sigGamPropInit;

aSig2 = sig2Params(1);
bSig2 = sig2Params(2);

if isempty(initDraws)
    rhoDraws = [muRho; NaN(burnIn+B-1,1)];
    gamDraws = [(aGam-1)/(aGam+bGam-2) ;NaN(burnIn+B-1,1)];
    sig2Draws = [bSig2/(aSig2+1); NaN(burnIn+B-1,1)];
else
    rhoDraws = [initDraws(1); NaN(burnIn+B-1,1)];
    gamDraws = [initDraws(2) ;NaN(burnIn+B-1,1)];
    sig2Draws = [initDraws(3); NaN(burnIn+B-1,1)];
end

F = [1 0; 0 rhoDraws(1)];
Q = [sqrt(gamDraws(1)*sig2Draws(1)) 0; 0 sqrt((1-gamDraws(1))*sig2Draws(1))];
H = [1 1];
Mdl = ssm(F,Q,H,'Mean0',[y(1);0],'Cov0',[1,-1;-1,1]);

statesSmoothedGuess = simsmooth(Mdl,y);
muDraws = [statesSmoothedGuess(:,1)'; NaN(burnIn+B-1,T)];
xDraws = [statesSmoothedGuess(:,2)'; NaN(burnIn+B-1,T)];

paramInd1 = logical([paramInd(1:end-1);0]);
paramInd2 = logical([0;paramInd(2:end)]);

while true

    for b = (1:burnIn-1)

        % Draw rho
        
        xTemp = xDraws(b,paramInd1);
        varRhoPosterior = 1/(1/varRho + sum((xTemp(keepInd)./sqrt((1-gamDraws(b))*sig2Draws(b))).^2));
        xTemp = xDraws(b,paramInd1).*xDraws(b,paramInd2);
        muRhoPosterior = varRhoPosterior*(muRho/varRho + sum(xTemp(keepInd))/...
            ((1-gamDraws(b))*sig2Draws(b)));
        rhoDraws(b+1) = normrnd(muRhoPosterior,sqrt(varRhoPosterior));

        % Draw gamma

        gamProp = normrnd(gamDraws(b),sigGamProp);
        xTemp = xDraws(b,paramInd2)-rhoDraws(b+1)*xDraws(b,paramInd1);
        muTemp = diff(muDraws(b,paramInd));

        logLCur = -(T1-1)*log(sig2Draws(b)) - (T1-1)/2*log(gamDraws(b)) - ...
            (T1-1)/2*log(1-gamDraws(b)) - sum((muTemp(keepInd)./sqrt(sig2Draws(b)*gamDraws(b))).^2)/2 - ...
            sum((xTemp(keepInd)./sqrt(sig2Draws(b)*(1-gamDraws(b)))).^2)/2 + ...
            logGamPrior(gamDraws(b));
        logLProp = -(T1-1)*log(sig2Draws(b)) - (T1-1)/2*log(gamProp) - ...
            (T1-1)/2*log(1-gamProp) - sum((muTemp(keepInd)./sqrt(sig2Draws(b)*gamProp)).^2)/2 - ...
            sum((xTemp(keepInd)./sqrt(sig2Draws(b)*(1-gamProp))).^2)/2 + ...
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

        aSig2Posterior = aSig2 + T1 - 1;
        bSig2Posterior = bSig2 + sum((muTemp(keepInd)./sqrt(gamDraws(b+1))).^2)/2 + ...
            sum((xTemp(keepInd)./sqrt(1-gamDraws(b+1))).^2)/2;
        sig2Draws(b+1) = 1/gamrnd(aSig2Posterior,1/bSig2Posterior);

        % Draw mu and x

        F = [1 0; 0 rhoDraws(b+1)];
        Q = [sqrt(gamDraws(b+1)*sig2Draws(b+1)) 0; 0 sqrt((1-gamDraws(b+1))*sig2Draws(b+1))];
        H = [1 1];
        Mdl = ssm(F,Q,H,'Mean0',[y(1);0],'Cov0',[1,-1;-1,1]);
        smoothedDraws = simsmooth(Mdl,y);
        muDraws(b+1,:) = smoothedDraws(:,1)';
        xDraws(b+1,:) = smoothedDraws(:,2)';

    end

    %mean(gamAcceptance)
    if mean(gamAcceptance) < 0.25
        sigGamProp = 0.8*sigGamProp;
        gamAcceptance = zeros(burnIn,1);
    elseif mean(gamAcceptance) > 0.4
        sigGamProp = 1.2*sigGamProp;
        gamAcceptance = zeros(burnIn,1);
    else
        break
    end

end

for b = (burnIn:burnIn+B-1)
    
    % Draw rho

    xTemp = xDraws(b,paramInd1);
        varRhoPosterior = 1/(1/varRho + sum((xTemp(keepInd)./sqrt((1-gamDraws(b))*sig2Draws(b))).^2));
    xTemp = xDraws(b,paramInd1).*xDraws(b,paramInd2);
    muRhoPosterior = varRhoPosterior*(muRho/varRho + sum(xTemp(keepInd))/...
        ((1-gamDraws(b))*sig2Draws(b)));
    rhoDraws(b+1) = normrnd(muRhoPosterior,sqrt(varRhoPosterior));

    % Draw gamma

    gamProp = normrnd(gamDraws(b),sigGamProp);
    xTemp = xDraws(b,paramInd2)-rhoDraws(b+1)*xDraws(b,paramInd1);
    muTemp = diff(muDraws(b,paramInd));

    logLCur = -(T1-1)*log(sig2Draws(b)) - (T1-1)/2*log(gamDraws(b)) - ...
        (T1-1)/2*log(1-gamDraws(b)) - sum((muTemp(keepInd)./sqrt(sig2Draws(b)*gamDraws(b))).^2)/2 - ...
        sum((xTemp(keepInd)./sqrt(sig2Draws(b)*(1-gamDraws(b)))).^2)/2 + ...
        logGamPrior(gamDraws(b));
    logLProp = -(T1-1)*log(sig2Draws(b)) - (T1-1)/2*log(gamProp) - ...
        (T1-1)/2*log(1-gamProp) - sum((muTemp(keepInd)./sqrt(sig2Draws(b)*gamProp)).^2)/2 - ...
        sum((xTemp(keepInd)./sqrt(sig2Draws(b)*(1-gamProp))).^2)/2 + ...
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

    aSig2Posterior = aSig2 + T1 - 1;
    bSig2Posterior = bSig2 + sum((muTemp(keepInd)./sqrt(gamDraws(b+1))).^2)/2 + ...
        sum((xTemp(keepInd)./sqrt(1-gamDraws(b+1))).^2)/2;
    sig2Draws(b+1) = 1/gamrnd(aSig2Posterior,1/bSig2Posterior);

    % Draw mu and x

    F = [1 0; 0 rhoDraws(b+1)];
    Q = [sqrt(gamDraws(b+1)*sig2Draws(b+1)) 0; 0 sqrt((1-gamDraws(b+1))*sig2Draws(b+1))];
    H = [1 1];
    Mdl = ssm(F,Q,H,'Mean0',[y(1);0],'Cov0',[1,-1;-1,1]);
    smoothedDraws = simsmooth(Mdl,y);
    muDraws(b+1,:) = smoothedDraws(:,1)';
    xDraws(b+1,:) = smoothedDraws(:,2)';

end

rhoDraws = rhoDraws(burnIn+1:end);
gamDraws = gamDraws(burnIn+1:end);
sig2Draws = sig2Draws(burnIn+1:end);
muDraws = muDraws(burnIn+1:end,:);
xDraws = xDraws(burnIn+1:end,:);

sigDraws = sqrt(sig2Draws);