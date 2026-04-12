% Epoch_Analysis_Plot.m
% Visualizes the Posterior Distributions of Gamma across historical epochs

clear; clc; close all;

disp('Loading Epoch Data...');

% Loads the epochs
load('gdpResults_AllCoefficients_Levels5.mat', 'gamDraws');
gam_Baseline = gamDraws;
load('gdpResults_GreatMod.mat', 'gamDraws');
gam_GreatMod = gamDraws;
load('gdpResults_ModernEra.mat', 'gamDraws');
gam_ModernEra = gamDraws;

disp('Generating Posterior Distribution Chart...');

% Calculates kernel density estimates 
% The (:) operator flattens any 3D arrays into a 1D vector so ksdensity doesn't crash
[f_base, xi_base] = ksdensity(gam_Baseline(:));
[f_gmod, xi_gmod] = ksdensity(gam_GreatMod(:));
[f_mod, xi_mod]  = ksdensity(gam_ModernEra(:));

% Creates the figure
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
hold on;

% Plots the distributions 
plot(xi_base, f_base, 'k-', 'LineWidth', 4, 'DisplayName', 'Baseline (1947 - 2019)');
plot(xi_gmod, f_gmod, 'b--', 'LineWidth', 4, 'DisplayName', 'Great Moderation (1984 - 2007)');
plot(xi_mod, f_mod, 'r:', 'LineWidth', 4, 'DisplayName', 'Modern Era (1984 - 2019)');
hold off;
title('Posterior Distribution of \gamma (Learning Parameter) Across Epochs', 'FontSize', 22, 'FontWeight', 'bold');
xlabel('Value of \gamma (Variance Share of Permanent Component)', 'FontSize', 18);
ylabel('Probability Density', 'FontSize', 18);
legend('Location', 'northeast', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 16);

% Saves the figure 
saveas(gcf, 'Epoch_Analysis_Gamma.png');
disp('Chart saved as Epoch_Analysis_Gamma.png!');