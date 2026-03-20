%% True model parameters
clear
clc

rho = 0.95;
gam = 0.3;
sig = 0.5;
trendInd = 0;

rho2=0.1;
gam2=0.01;

%% Prior parameters (modify to match existing results)

rhoParams_unb = [0.95; 0.1^2];% mean and variance of Normal prior for rho
rhoParams_dwb= [0.4,0.1^2];
rhoParams_uwb= [0.95,0.1^2];

%gamMoments = [0.01; 0.01]; % mode and standard deviation of Beta prior for gamma
%gamFun = @(x) [(x(1)-1)/(sum(x)-2); sqrt(prod(x)/((sum(x)^2)*(sum(x)+1)))]; % moments as function of Beta distribution parameters
%gamParams = fminunc(@(x) 1000*sum((gamFun(x)-gamMoments).^2),[2;2]); % solve for Beta distribution parameters associated with moments
gamParams_unb = [9.052;19.788];
gamParams_dwb = [2.34, 26.5];
gamParams_uwb= [9.052, 19.788];


sig2Params = [1.25; 0.5625]; % Inverse Gamma distribution prior parameters for conditional variance, do not change

invGamPdf = @(x,a,b) ((b^a)/gamma(a))*(x.^(-a-1)).*exp(-b./x).*(x>0);

% Distributions
rhoPrior_unb = @(x) normpdf(x,rhoParams_unb(1),sqrt(rhoParams_unb(2)));
gamPrior_unb = @(x) betapdf(x,gamParams_unb(1),gamParams_unb(2));

rhoPrior_dwb = @(x) normpdf(x,rhoParams_dwb(1),sqrt(rhoParams_dwb(2)));
gamPrior_dwb = @(x) betapdf(x,gamParams_dwb(1),gamParams_dwb(2));

rhoPrior_uwb = @(x) normpdf(x,rhoParams_uwb(1),sqrt(rhoParams_uwb(2)));
gamPrior_uwb = @(x) betapdf(x,gamParams_uwb(1),gamParams_uwb(2));

sigPrior = @(x) invGamPdf(x,sig2Params(1),sig2Params(2));

priorPlot1 = linspace(0,1,1000)';
priorPlot2 = linspace(0,2,1000)';
figure
subplot(3,3,1)
plot(priorPlot1,rhoPrior_unb(priorPlot1),'k-','LineWidth',3)
xline(rho,'Color','#808080','LineWidth',3);
xt = xticks;
title('\rho')
set(gca,'FontSize',20)
axis tight

subplot(3,3,2)
plot(priorPlot1,gamPrior_unb(priorPlot1),'k-','LineWidth',3)
xline(gam,'Color','#808080','LineWidth',3);
title('\gamma')
set(gca,'FontSize',20)
axis tight

subplot(3,3,3)
plot(priorPlot2,sigPrior(priorPlot2),'k-','LineWidth',3)
xline(sig^2,'Color','#808080','LineWidth',3);
title('\sigma^2')
set(gca,'FontSize',20)
axis tight


subplot(3,3,4)
plot(priorPlot1,rhoPrior_dwb(priorPlot1),'k-','LineWidth',3)
xline(rho,'Color','#808080','LineWidth',3);
xt = xticks;
title('\rho')
set(gca,'FontSize',20)
axis tight

subplot(3,3,5)
plot(priorPlot1,gamPrior_dwb(priorPlot1),'k-','LineWidth',3)
xline(gam,'Color','#808080','LineWidth',3);
title('\gamma')
set(gca,'FontSize',20)
axis tight

subplot(3,3,6)
plot(priorPlot2,sigPrior(priorPlot2),'k-','LineWidth',3)
xline(sig^2,'Color','#808080','LineWidth',3);
title('\sigma^2')
set(gca,'FontSize',20)
axis tight

subplot(3,3,7)
plot(priorPlot1,rhoPrior_uwb(priorPlot1),'k-','LineWidth',3)
xline(rho2,'Color','#808080','LineWidth',3);
xt = xticks;
title('\rho')
set(gca,'FontSize',20)
axis tight

subplot(3,3,8)
plot(priorPlot1,gamPrior_uwb(priorPlot1),'k-','LineWidth',3)
xline(gam2,'Color','#808080','LineWidth',3);
title('\gamma')
set(gca,'FontSize',20)
axis tight

subplot(3,3,9)
plot(priorPlot2,sigPrior(priorPlot2),'k-','LineWidth',3)
xline(sig^2,'Color','#808080','LineWidth',3);
title('\sigma^2')
set(gca,'FontSize',20)
axis tight
