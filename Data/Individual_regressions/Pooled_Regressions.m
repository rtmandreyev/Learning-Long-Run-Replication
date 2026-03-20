clear
clc
%% Bias regressions:
% load data 
% In these tables, column 3 is the dependent variable of the bias
% regressions; column 1 is the forecaster ID, column 2 is the quarter
% number: for example: 1981Q3 is number 1 and 2018 Q4 is number 150

Bias_h1data=readtable("Bias_data\Bias_data_mat_h1.csv");
Bias_h2data=readtable("Bias_data\Bias_data_mat_h2.csv");
Bias_h3data=readtable("Bias_data\Bias_data_mat_h3.csv");
Bias_h4data=readtable("Bias_data\Bias_data_mat_h4.csv");



%% Load Autocorrelated errors data:
% In these tables, column 1 is the forecaster ID, column 2 is the quarter
% number, column 3 is the dependent variable, column 4 is the independent
% variable of the autocorrelated errors regressions

AR1Forecasterdata=readtable("AR_data\AR1_Forecaster_data.csv");
AR2Forecasterdata=readtable("AR_data\AR2_Forecaster_data.csv");
AR3Forecasterdata=readtable("AR_data\AR3_Forecaster_data.csv");
AR4Forecasterdata=readtable("AR_data\AR4_Forecaster_data.csv");

%% Load Mincer Zarnowitz data:
% In these tables, column 1 is the forecaster ID, column 2 is the quarter
% number, column 4 is the dependent variable, column 3 is the independent
% variable of the Mincer Zarnowitz regressions

MZ_h1data=readtable("MZ_data\MZ_data_mat_h1.csv");
MZ_h2data=readtable("MZ_data\MZ_data_mat_h2.csv");
MZ_h3data=readtable("MZ_data\MZ_data_mat_h3.csv");
MZ_h4data=readtable("MZ_data\MZ_data_mat_h4.csv");

%% Load Coibion Gorodnichenko data:
CG_h1data=readtable("CG_data\CG_h1_data.csv");
CG_h2data=readtable("CG_data\CG_h2_data.csv");
CG_h3data=readtable("CG_data\CG_h3_data.csv");

%% Compute forecast anomalies
maxHorizon=4;
biasMatModel = NaN(maxHorizon+1,1);
biasMatModel_se = biasMatModel;
biasMatModel_p = biasMatModel;

arMatModel = NaN(maxHorizon+1,2);
arMatModel_se = arMatModel;
arMatModel_p = arMatModel;

mzMatModel = NaN(maxHorizon+1,2);
mzMatModel_se = mzMatModel;
mzMatModel_p = mzMatModel;

cgMatModel = NaN(maxHorizon+1,2);
cgMatModel_se = cgMatModel;
cgMatModel_p = cgMatModel;

%% Bias Regressions
temp=regstats2(table2array(Bias_h1data(:,3)),ones(size(Bias_h1data,1),1),...
        'onlydata',{'beta','hac'});
biasMatModel(1)=temp.beta
biasMatModel_se(1)=temp.hac.se
biasMatModel_p(1)=temp.hac.pval

temp=regstats2(table2array(Bias_h2data(:,3)),ones(size(Bias_h2data,1),1),...
        'onlydata',{'beta','hac'});
biasMatModel(2)=temp.beta;
biasMatModel_se(2)=temp.hac.se;
biasMatModel_p(2)=temp.hac.pval;

temp=regstats2(table2array(Bias_h3data(:,3)),ones(size(Bias_h3data,1),1),...
        'onlydata',{'beta','hac'});
biasMatModel(3)=temp.beta;
biasMatModel_se(3)=temp.hac.se;
biasMatModel_p(3)=temp.hac.pval;

temp=regstats2(table2array(Bias_h4data(:,3)),ones(size(Bias_h4data,1),1),...
        'onlydata',{'beta','hac'});
biasMatModel(4)=temp.beta;
biasMatModel_se(4)=temp.hac.se;
biasMatModel_p(4)=temp.hac.pval;








%% Autocorrelated errors regressions

temp=regstats2(table2array(AR1Forecasterdata(:,3)),table2array(AR1Forecasterdata(:,4)), 'linear',{'beta', 'hac'})
arMatModel(1,:) = temp.beta;
arMatModel_se(1,:) = temp.hac.se;
arMatModel_p(1,:) = temp.hac.pval;
ar1_results_p=[temp.beta(2) temp.hac.se(2) temp.hac.pval(2)];

temp=regstats2(table2array(AR2Forecasterdata(:,3)),table2array(AR2Forecasterdata(:,4)), 'linear',{'beta', 'hac'})
arMatModel(2,:) = temp.beta;
arMatModel_se(2,:) = temp.hac.se;
arMatModel_p(2,:) = temp.hac.pval;
ar2_results_p=[temp.beta(2) temp.hac.se(2) temp.hac.pval(2)];

temp=regstats2(table2array(AR3Forecasterdata(:,3)),table2array(AR3Forecasterdata(:,4)), 'linear',{'beta', 'hac'})
arMatModel(3,:) = temp.beta;
arMatModel_se(3,:) = temp.hac.se;
arMatModel_p(3,:) = temp.hac.pval;
ar3_results_p=[temp.beta(2) temp.hac.se(2) temp.hac.pval(2)];

temp=regstats2(table2array(AR4Forecasterdata(:,3)),table2array(AR4Forecasterdata(:,4)), 'linear',{'beta', 'hac'})
arMatModel(4,:) = temp.beta;
arMatModel_se(4,:) = temp.hac.se;
arMatModel_p(4,:) = temp.hac.pval;
ar4_results_p=[temp.beta(2) temp.hac.se(2) temp.hac.pval(2)];


%% Mincer Zarnowitz regressions
 temp=regstats2(table2array(MZ_h1data(:,4)),table2array(MZ_h1data(:,3)),...
        'linear',{'beta','hac'});
 mzMatModel(1,:) = temp.beta;
 mzMatModel_se(1,:) = temp.hac.se;
 mzMatModel_p(1,:) = [temp.hac.pval(1),2*normcdf(-abs((temp.beta(2)-1)/temp.hac.se(2)))];


 temp=regstats2(table2array(MZ_h2data(:,4)),table2array(MZ_h2data(:,3)),...
        'linear',{'beta','hac'});
 mzMatModel(2,:) = temp.beta;
 mzMatModel_se(2,:) = temp.hac.se;
 mzMatModel_p(2,:) = [temp.hac.pval(1),2*normcdf(-abs((temp.beta(2)-1)/temp.hac.se(2)))];

 temp=regstats2(table2array(MZ_h3data(:,4)),table2array(MZ_h3data(:,3)),...
        'linear',{'beta','hac'});
 mzMatModel(3,:) = temp.beta;
 mzMatModel_se(3,:) = temp.hac.se;
 mzMatModel_p(3,:) = [temp.hac.pval(1),2*normcdf(-abs((temp.beta(2)-1)/temp.hac.se(2)))];

 temp=regstats2(table2array(MZ_h4data(:,4)),table2array(MZ_h4data(:,3)),...
        'linear',{'beta','hac'});
 mzMatModel(4,:) = temp.beta;
 mzMatModel_se(4,:) = temp.hac.se;
 mzMatModel_p(4,:) = [temp.hac.pval(1),2*normcdf(-abs((temp.beta(2)-1)/temp.hac.se(2)))];







%% Coibion Gorodnichenko regressions

temp=regstats2(table2array(CG_h1data(:,3)),table2array(CG_h1data(:,4)), 'linear',{'beta','hac'});

cgMatModel(1,:)=temp.beta;
cgMatModel_se(1,:)=temp.hac.se;
cgMatModel_p(1,:)=temp.hac.pval;


temp=regstats2(table2array(CG_h2data(:,3)),table2array(CG_h2data(:,4)), 'linear',{'beta','hac'});

cgMatModel(2,:)=temp.beta;
cgMatModel_se(2,:)=temp.hac.se;
cgMatModel_p(2,:)=temp.hac.pval;


temp=regstats2(table2array(CG_h3data(:,3)),table2array(CG_h3data(:,4)), 'linear',{'beta','hac'});

cgMatModel(3,:)=temp.beta;
cgMatModel_se(3,:)=temp.hac.se;
cgMatModel_p(3,:)=temp.hac.pval;



%%
Bias_results = [biasMatModel(1:4) biasMatModel_se(1:4) biasMatModel_p(1:4)]'
AR_results= [arMatModel(1:4,2) arMatModel_se(1:4,2) arMatModel_p(1:4,2)]'
MZ_results= [mzMatModel(1:4,2) mzMatModel_se(1:4,2) mzMatModel_p(1:4,2)]'  
CG_results=[cgMatModel(1:3,2) cgMatModel_se(1:3,2) cgMatModel_p(1:3,2)]'