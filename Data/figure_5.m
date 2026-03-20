%% Import data from spreadsheet


%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:B295";

% Specify column names and types
opts.VariableNames = ["Year", "ConsoleRate"];
opts.VariableTypes = ["double", "double"];

% Import the data
UKConsoleRate = readtable("C:\Users\shara\OneDrive\Desktop\Replication_Package_Learning\Data\UKConsoleRate.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts
%%
UK_Consol_data= table2array(UKConsoleRate);
%%
fig5 = figure(5);

plot(UK_Consol_data(:,1),UK_Consol_data(:,2),'-k','LineWidth',3)
xline(1951,'--')
set(findall(gcf,'-property','FontSize'),'FontSize',30)
fig5.WindowState = 'maximized';
pause(1)
hold on
