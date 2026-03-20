% Single-click wrapper for Farmer, Nakamura, and Steinsson (2024)

clear; clc; close all;

% =========================================================================
% MASTER CONTROL PANEL
% INSTRUCTOR NOTE: Set to 'true' to run the hours-long estimations.
% =========================================================================
run_full_tbill_estimation = true; % Set false for subsequent runs
run_full_gdp_estimation   = true; % Set false for subsequent runs
run_full_mc_estimation    = false; % Pre-computed files provided by authors

save('system_toggles.mat', 'run_full_tbill_estimation', 'run_full_gdp_estimation', 'run_full_mc_estimation');
% =========================================================================

disp('======================================================');
disp('STARTING REPLICATION: Learning About the Long Run');
disp('======================================================');

if exist('Tables_Output.txt', 'file'), delete('Tables_Output.txt'); end
diary('Tables_Output.txt');

% --- FOLDER 1: DATA ---
disp('Navigating to Data folder...');
cd('Data');
data_regs; data_figs;

copyfile('figure_5.m', 'figure_5_backup.m');
txt = fileread('figure_5.m');
txt = strrep(txt, '"C:\Users\shara\OneDrive\Desktop\Replication_Package_Learning\Data\UKConsoleRate.xlsx"', '"UKConsoleRate.xlsx"');
fid = fopen('figure_5.m', 'w'); fwrite(fid, txt, '*char'); fclose(fid);

figure_5;
copyfile('figure_5_backup.m', 'figure_5.m'); delete('figure_5_backup.m');
cd('..'); exportOpenFigures(pwd, 'Data_Fig');

% --- FOLDER 2: T-BILL RESULTS ---
disp('Navigating to T-bill results folder...');
cd('T-bill results');
load('../system_toggles.mat'); 

if run_full_tbill_estimation
    disp('Running tblPriorEstimationSmooth.m (Full Estimation - Hours)...');
    tblPriorEstimationSmooth;
else
    disp('Loading pre-computed 2GB file...');
end

copyfile('tblRegressions.m', 'tblRegressions_backup.m');
txt = fileread('tblRegressions.m');
txt = strrep(txt, 'tblResults_LookAhead_SmoothL.mat', 'tblResults_SlopeCoefficients_Level_SmoothL.mat');
fid = fopen('tblRegressions.m', 'w'); fwrite(fid, txt, '*char'); fclose(fid);

copyfile('tblFigs.m', 'tblFigs_backup.m');
txt = fileread('tblFigs.m');
txt = strrep(txt, 'tblResults_LookAhead_SmoothL.mat', 'tblResults_SlopeCoefficients_Level_SmoothL.mat');
txt = strrep(txt, 'gam1Params', 'gamParams'); txt = strrep(txt, 'rho1Params', 'rhoParams'); 
txt = strrep(txt, 'gam2Params', 'gamParams'); txt = strrep(txt, 'rho2Params', 'rhoParams'); 
fid = fopen('tblFigs.m', 'w'); fwrite(fid, txt, '*char'); fclose(fid);

tblRegressions; tblFigs;
copyfile('tblRegressions_backup.m', 'tblRegressions.m'); delete('tblRegressions_backup.m');
copyfile('tblFigs_backup.m', 'tblFigs.m'); delete('tblFigs_backup.m');

cd('..'); exportOpenFigures(pwd, 'Tbill_Fig');
disp('Skipping tblPriorEstimation_Break.m (Figure 9) due to 18.4GB RAM limitation...');

% --- FOLDER 3: GDP RESULTS ---
disp('Navigating to GDP Results folder...');
cd('GDP Results');
load('../system_toggles.mat'); 

if run_full_gdp_estimation
    disp('Running gdpPriorEstimation.m (Full Estimation - Hours)...');
    gdpPriorEstimation;
else
    disp('Loading pre-computed file...');
end

copyfile('gdpFigs.m', 'gdpFigs_backup.m');
txt = fileread('gdpFigs.m');
load_anchor = 'load(''Gdp_Regression_Stats.mat'')';
patch = [newline 'try, if isdatetime(dateVec), dateVec = year(dateVec) + (month(dateVec)-1)/12; elseif iscell(dateVec)||isstring(dateVec), tmp=datetime(dateVec); dateVec=year(tmp)+(month(tmp)-1)/12; end, catch, end' newline];
txt = strrep(txt, load_anchor, [load_anchor patch]);
fid = fopen('gdpFigs.m', 'w'); fwrite(fid, txt, '*char'); fclose(fid);

gdpRegressions; gdpFigs;
copyfile('gdpFigs_backup.m', 'gdpFigs.m'); delete('gdpFigs_backup.m');

cd('..'); exportOpenFigures(pwd, 'GDP_Fig');

% --- FOLDER 4: MONTE CARLO ---
disp('Navigating to Monte_Carlo folder...');
cd('Monte_Carlo');
load('../system_toggles.mat'); 

if run_full_mc_estimation
    ucSimulation2;
else
    disp('Skipping ucSimulation2.m. Loading pre-computed files...');
end

% Backup the MC script 
copyfile('MonteCarlo_Results.m', 'MonteCarlo_Results_backup.m');

% Run 1: Downward Biased
disp('Generating Figure 14 and Tables 6 & 7 (Downward Biased)...');
txt = fileread('MonteCarlo_Results_backup.m'); 
fid = fopen('MonteCarlo_Results.m', 'w'); fwrite(fid, txt, '*char'); fclose(fid);
MonteCarlo_Results;
cd('..'); exportOpenFigures(pwd, 'MC_Fig14'); cd('Monte_Carlo');

% Run 2: Unbiased
disp('Generating Tables 6 & 7 (Unbiased)...');
txt = fileread('MonteCarlo_Results_backup.m'); 
txt = strrep(txt, 'data = ''New_DownwardBiasedPriors'';', 'data = ''New_UnBiasedPriors'';');
fid = fopen('MonteCarlo_Results.m', 'w'); fwrite(fid, txt, '*char'); fclose(fid);
MonteCarlo_Results; close all;

% Run 3: Upward Biased
disp('Generating Tables 6 & 7 (Upward Biased)...');
txt = fileread('MonteCarlo_Results_backup.m'); 
txt = strrep(txt, 'data = ''New_DownwardBiasedPriors'';', 'data = ''New_UpwardBiasedPriors'';');
fid = fopen('MonteCarlo_Results.m', 'w'); fwrite(fid, txt, '*char'); fclose(fid);
MonteCarlo_Results; close all;

% Restore original file
copyfile('MonteCarlo_Results_backup.m', 'MonteCarlo_Results.m'); delete('MonteCarlo_Results_backup.m');

% Final Figures
disp('Generating Figures 13 and 15...');
Figure_13; cd('..'); exportOpenFigures(pwd, 'MC_Fig13'); cd('Monte_Carlo');
Figure_15; cd('..'); exportOpenFigures(pwd, 'MC_Fig15');

% Cleanup
delete('system_toggles.mat');
diary off;
disp('======================================================');
disp('REPLICATION COMPLETE. All Figures and Tables Generated.');
disp('======================================================');

function exportOpenFigures(saveDir, prefix)
    figs = findobj('Type', 'figure');
    for i = 1:length(figs)
        figHandle = figs(i);
        if isvalid(figHandle)
            fileName = fullfile(saveDir, sprintf('%s_%d.png', prefix, figHandle.Number));
            print(figHandle, fileName, '-dpng', '-r300');
            close(figHandle);
        end
    end
    drawnow;
end