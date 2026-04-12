% gdp_ModernEra.m
% Epoch Analysis: The Modern Era (1984 - 2019)

clear; clc;
load('Gdp_Data.mat');

disp('Configuring Data for the Modern Era (1984 - 2019)...');

% Finds the dates for the epoch
epoch_idx = year(dateVec) >= 1984 & year(dateVec) <= 2019;

% Extracts only the final, fully revised data column
y_final_vintage = gdpData(epoch_idx, end);

% Removes NaNs and creates the clean matrix
y_final_vintage = y_final_vintage(~isnan(y_final_vintage));
T_epoch = length(y_final_vintage);
gdpData_clean = repmat(y_final_vintage, 1, T_epoch + 10); 

% Adjusts the MCMC starting markers
estimationStart = T_epoch; 
vintageStart = 1; 

% Authors' original parameters
burnIn = 25000;
B = 50000;
stepSize = 4;
seed = 15284;

thetOpt = [0.8888; 0.0028; -0.9129; 0.0025; 45.9721; 39.3704];
rhoParams = thetOpt(1:2);
rho2Params = thetOpt(3:4);
gamParams = thetOpt(5:6);
muParams = [0.01; 0.01^2];
sig2Moments = [0.015^2; 0.01^2];
sig2Params = [sig2Moments(1)^2/(sig2Moments(2)^2)+2; NaN];
sig2Params(2) = sig2Moments(1)*(sig2Params(1)-1);

% Runs the MCM engine
disp(['Starting MCMC Estimation for ', num2str(T_epoch), ' quarters...']);
[muDraws,gamDraws,rho1Draws,rho2Draws,sigDraws,~,~,~,~,~,~] = ...
    recursiveGibbsUC_GDP(log(gdpData_clean),muParams,gamParams,rhoParams,...
    rho2Params,sig2Params,burnIn,B,estimationStart,vintageStart,...
    stepSize,true(T_epoch,1),seed,1);

% Saves the required results
disp('Saving Results...');
save('gdpResults_ModernEra.mat', 'gamDraws', 'muDraws', 'rho1Draws', 'rho2Draws', 'sigDraws');
disp('Modern Era Analysis Complete!');