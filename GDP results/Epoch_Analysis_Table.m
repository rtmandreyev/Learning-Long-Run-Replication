% Epoch_Analysis_Table.m
% Generates a statistical summary table for the epoch analysis

clear; clc;

disp('Loading Epoch Data...');

% Loads the epochs
load('gdpResults_AllCoefficients_Levels5.mat', 'gamDraws');
gam_Baseline = gamDraws(:);
load('gdpResults_GreatMod.mat', 'gamDraws');
gam_GreatMod = gamDraws(:);
load('gdpResults_ModernEra.mat', 'gamDraws');
gam_ModernEra = gamDraws(:);

disp('Calculating Summary Statistics...');

% Helper function with the 'omitnan' flag to prevent crashes
calcStats = @(x) [mean(x, 'omitnan'), median(x, 'omitnan'), std(x, 'omitnan'), prctile(x, 5), prctile(x, 95)];

stats_Base = calcStats(gam_Baseline);
stats_GMod = calcStats(gam_GreatMod);
stats_Mod = calcStats(gam_ModernEra);

% Constructs the Table
Epochs = {'Baseline (1947-2019)'; 'Great Moderation (1984-2007)'; 'Modern Era (1984-2019)'};
Mean = [stats_Base(1); stats_GMod(1); stats_Mod(1)];
Median = [stats_Base(2); stats_GMod(2); stats_Mod(2)];
Std_Dev = [stats_Base(3); stats_GMod(3); stats_Mod(3)];
Pct_5th = [stats_Base(4); stats_GMod(4); stats_Mod(4)];
Pct_95th = [stats_Base(5); stats_GMod(5); stats_Mod(5)];
SummaryTable = table(Mean, Median, Std_Dev, Pct_5th, Pct_95th, 'RowNames', Epochs);

disp(' ');
disp('================================================================================');
disp('      Table 1: Posterior Estimates of Gamma Across Historical Epochs');
disp('================================================================================');
disp(SummaryTable);
disp('================================================================================');