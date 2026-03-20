clear,clc;

tic;
savename = 'simulationGibbs_New_DownwardBiasedPriors.mat';

%% True model parameters

rho = 0.95;
gam = 0.3;
sig = 0.5;

'Prior parameters (modify to match existing results)'

rhoParams = [0.4; 0.1^2]; % mean and variance of Normal prior for rho
gamMoments = [0.01; 0.05]; % mode and standard deviation of Beta prior for gamma
gamFun = @(x) [(x(1)-1)/(sum(x)-2); sqrt(prod(x)/((sum(x)^2)*(sum(x)+1)))]; % moments as function of Beta distribution parameters
gamParams = fminunc(@(x) 1000*sum((gamFun(x)-gamMoments).^2),[2;2]); % solve for Beta distribution parameters associated with moments
gamParams = [2.34;26.5];
sig2Params = [1.25; 0.5625]; % Inverse Gamma distribution prior parameters for conditional variance, do not change


'Set simulation parameters and initialize structures'

load('Tbl_Data.mat','tblData','estimationStart','T','spfForecastsQ','yieldData')

% number of bootstrap replications

%Ethan's Mac
%poolobj = parpool('local',2);
%B = 2;

%Cluster
cl = parcluster();
cl.NumWorkers = str2num(getenv('SLURM_CPUS_PER_TASK'));
poolobj = cl.parpool(str2num(getenv('SLURM_CPUS_PER_TASK')));
B = 500; %Number of simulations

xSamples = [sig/sqrt(1-rho^2)*ones(1,B); NaN(T-1,B)];
muSamples = [tblData(1)-xSamples(1,:); NaN(T-1,B)];

'Simulate bootstrap samples'

uSim = sqrt(gam)*sig*randn(T,B); % shocks to mu
vSim = sqrt(1-gam)*sig*randn(T,B); % shocks to x

for t = 1:(T-1)

    muSamples(t+1,:) = muSamples(t,:) + uSim(t+1,:);
    xSamples(t+1,:) = rho*xSamples(t,:) + vSim(t+1,:);

end

tblSamples = muSamples + xSamples;

'Main bootstrap simulations'

% These parameters are set to 25,000 for our data excercise but I don't
% think it would be feasible to run repeated simulations with those numbers
% Tentatively set to 5,000 each

burnIn = 5000; % number of draws to discard as burn-in in Gibbs Sampler
numKeep = 5000; % number of draws to keep after burn-in period in Gibbs Sampler
stepSize = 4; % re-estimate parameters every 4 periods, can also play with this if simulations are too costly
maxYieldHorizon = 40;

% Mean parameter estimates

rhoMeanSim = NaN(T,B);
gamMeanSim = NaN(T,B);
sigMeanSim = NaN(T,B);
meanForecastsSim = NaN(T,maxYieldHorizon,B);
parfor b = 1:B
	
	strcat('Simulation:',num2str(b),' Start')	
    % Gibbs Sampler

    [rhoDraws,gamDraws,sigDraws,~,~,~,~,yForecasts] = ...
        recursiveGibbsUC(tblSamples(:,b),rhoParams,gamParams,sig2Params,...
        burnIn,numKeep,estimationStart,stepSize,ones(T,1),'shuffle');

    % Save parameter estiamtes

    rhoMeanSim(:,b) = nanmean(rhoDraws,2);
    gamMeanSim(:,b) = nanmean(gamDraws,2);
    sigMeanSim(:,b) = nanmean(sigDraws,2);

    % Align timing of forecasts and set sample same as SPF

    meanForecastsSim(:,:,b) = squeeze(nanmean(yForecasts,3));

    strcat('Simulation:',num2str(b),' End')
end
delete(poolobj);

save(savename,'rho','gam','sig','T','B','rhoParams','gamParams','sig2Params',...
    'rhoMeanSim','gamMeanSim','sigMeanSim','meanForecastsSim','burnIn',...
    'numKeep','stepSize','maxYieldHorizon','muSamples','xSamples');
toc;
