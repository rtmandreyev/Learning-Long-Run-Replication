clear
clc
%%
Final_Individual_For = readtable("Individual_data_n.csv");
data_forecast_match = readtable("data_IFC_Agg_n.csv");
CGdataID=readtable("CG_data_ID.csv");
%% Create the data required to run the regeressions
data_forecast_agg=table2array(data_forecast_match);
Individual_data=table2array(Final_Individual_For(:,1:8));

Individual_data= [Individual_data table2array(Final_Individual_For(:,10))];
%%
size(data_forecast_agg,1);

unique_ID=unique(Individual_data(:,3)); %correct

bias_coefs=zeros(size(unique_ID,1),4);
mz_coefs_alpha=bias_coefs;
mz_coefs_beta=bias_coefs;

Bias_data_mat=zeros(size(data_forecast_agg,1),5);
Bias_data_mat(:,1)=Individual_data(:,3);

MZ_aggregate_data_mat=zeros(size(data_forecast_agg,1),5);
MZ_aggregate_data_mat(:,1)=Individual_data(:,3);
MZ_forecaster_data_mat=zeros(size(data_forecast_agg,1),5);
MZ_forecaster_data_mat(:,1)=Individual_data(:,3);
for i=2:5
    Bias_data_mat(:,i) = data_forecast_agg(:,i) - Individual_data(:,3+i);
    MZ_aggregate_data_mat(:,i)= data_forecast_agg(:,i);
    MZ_forecaster_data_mat(:,i)=Individual_data(:,3+i);
end

%%
data_forecast_merged=[Individual_data data_forecast_agg];
% Variables are: year ; quarter ; ID ; Nowcast ; 1Q Ahead ; 2Q Ahead; 3Q Ahead; 4Q Ahead; NowCD ; 1QAhead_CD ; 2Q_CD ; 3Q_CD ; 4Q_CD  

unique_ID(1)

ID_ind=find(Individual_data(:,3) == 2);

Bias_data_ind=Bias_data_mat(ID_ind,:);

MZ_aggdata_ind=MZ_aggregate_data_mat(ID_ind,:);
MZ_fcdata_ind= MZ_forecaster_data_mat(ID_ind,:);

%% Forecaster by forecaster bias and mincer zarnowitz regressions:

for j=1:4
    for i = 1:size(bias_coefs,1)
        ID_ind=find(Bias_data_mat(:,1)==unique_ID(i));
        Bias_data_ind=Bias_data_mat(ID_ind,:);
        MZ_aggdata_ind=MZ_aggregate_data_mat(ID_ind,:);
        MZ_fcdata_ind= MZ_forecaster_data_mat(ID_ind,:);
        if sum(~isnan(Bias_data_ind(:,j+1))) >=10
            temp=regstats2(Bias_data_ind(:,j+1),ones(size(Bias_data_ind,1),1),'onlydata',{'beta','hac'});
            bias_coefs(i,j)=temp.beta;
        end

        if sum(~isnan(MZ_fcdata_ind(:,j+1)))  >=10
            temp1=regstats2(MZ_aggdata_ind(:,j+1),MZ_fcdata_ind(:,j+1),'linear',{'beta','hac'});
            mz_coefs_alpha(i,j)=temp1.beta(1);
            mz_coefs_beta(i,j)=temp1.beta(2);
        end

    end
end
%% correct bias coefficients
bias_coefs_h1=bias_coefs(bias_coefs(:,1)~=0,1);


bias_coefs_h2=bias_coefs(bias_coefs(:,2)~=0,2);


bias_coefs_h3=bias_coefs(bias_coefs(:,3)~=0,3);

bias_coefs_h4=bias_coefs(bias_coefs(:,4)~=0,4);




%%
mz_coefbeta_g0_h1=mz_coefs_beta(mz_coefs_beta(:,1)~=0,1);
mz_coefbeta_g0_h2=mz_coefs_beta(mz_coefs_beta(:,2)~=0,2);
mz_coefbeta_g0_h3=mz_coefs_beta(mz_coefs_beta(:,3)~=0,3);
mz_coefbeta_g0_h4=mz_coefs_beta(mz_coefs_beta(:,4)~=0,4);








%% CG regressions

unique_ID_CG=table2array(unique(CGdataID(:,7)));
CG_Data_ID = table2array(CGdataID);
unique_index_CG=find(CG_Data_ID(:,7)==unique_ID_CG(4));

CG_Data_individual=CG_Data_ID(unique_index_CG,:);

CG_coefficients=zeros(length(unique_ID_CG),3);

size(CG_Data_individual,1);

%%
for i=1:length(unique_ID_CG)
    unique_index_CG=find(CG_Data_ID(:,7)==unique_ID_CG(i));
    CG_Data_individual=CG_Data_ID(unique_index_CG,:);
    if sum(~isnan(CG_Data_individual(:,1)))>9 && sum(~isnan(CG_Data_individual(:,2)))>9
        temp1=regstats2(CG_Data_individual(:,1),CG_Data_individual(:,2),'linear',{'beta','hac'});
        CG_coefficients(i,1)=temp1.beta(2);
    end
    
    if sum(~isnan(CG_Data_individual(:,3)))>9 && sum(~isnan(CG_Data_individual(:,4)))>9
        temp2=regstats2(CG_Data_individual(:,3),CG_Data_individual(:,4),'linear',{'beta','hac'});
        CG_coefficients(i,2)=temp2.beta(2);
    end
    
    if sum(~isnan(CG_Data_individual(:,5)))>9 && sum(~isnan(CG_Data_individual(:,6)))>9
        temp3=regstats2(CG_Data_individual(:,5),CG_Data_individual(:,6),'linear',{'beta','hac'});
        CG_coefficients(i,3)=temp3.beta(2);
    end
    


end
%%
cg_h1_coefs=CG_coefficients(CG_coefficients(:,1)~=0,1);
cg_h2_coefs=CG_coefficients(CG_coefficients(:,2)~=0,2);
cg_h3_coefs=CG_coefficients(CG_coefficients(:,3)~=0,3);

%%
AR1Forecasterdata=table2array(readtable("AR_data\AR1_Forecaster_data.csv"));
AR2Forecasterdata=table2array(readtable("AR_data\AR2_Forecaster_data.csv"));
AR3Forecasterdata=table2array(readtable("AR_data\AR3_Forecaster_data.csv"));
AR4Forecasterdata=table2array(readtable("AR_data\AR4_Forecaster_data.csv"));
%%
unique_ID_ar1=unique(AR1Forecasterdata(:,1));
ar1_coefficients=zeros(length(unique_ID_ar1),1);

for i=1:length(unique_ID_ar1)
    index_ind_ar1=find(AR1Forecasterdata(:,1)==unique_ID_ar1(i));
    ar1_fcdata_individual=AR1Forecasterdata(index_ind_ar1,:);
    temp=regstats2(ar1_fcdata_individual(:,3),ar1_fcdata_individual(:,4), 'linear',{'beta', 'hac'});
    ar1_coefficients(i)=temp.beta(2);
end

unique_ID_ar2=unique(AR2Forecasterdata(:,1));
ar2_coefficients=zeros(length(unique_ID_ar2),1);

for i=1:length(unique_ID_ar2)
    index_ind_ar2=find(AR2Forecasterdata(:,1)==unique_ID_ar2(i));
    ar2_fcdata_individual=AR2Forecasterdata(index_ind_ar2,:);
    temp=regstats2(ar2_fcdata_individual(:,3),ar2_fcdata_individual(:,4), 'linear',{'beta', 'hac'});
    ar2_coefficients(i)=temp.beta(2);
end

unique_ID_ar3=unique(AR3Forecasterdata(:,1));
ar3_coefficients=zeros(length(unique_ID_ar3),1);
for i=1:length(unique_ID_ar3)
    index_ind_ar3=find(AR3Forecasterdata(:,1)==unique_ID_ar3(i));
    ar3_fcdata_individual=AR3Forecasterdata(index_ind_ar3,:);
    temp=regstats2(ar3_fcdata_individual(:,3),ar3_fcdata_individual(:,4), 'linear',{'beta', 'hac'});
    ar3_coefficients(i)=temp.beta(2);
end

unique_ID_ar4=unique(AR4Forecasterdata(:,1));
ar4_coefficients=zeros(length(unique_ID_ar4),1);
for i=1:length(unique_ID_ar4)
    index_ind_ar4=find(AR4Forecasterdata(:,1)==unique_ID_ar4(i));
    ar4_fcdata_individual=AR4Forecasterdata(index_ind_ar4,:);
    temp=regstats2(ar4_fcdata_individual(:,3),ar4_fcdata_individual(:,4), 'linear',{'beta', 'hac'});
    ar4_coefficients(i)=temp.beta(2);
end
%%
[median(bias_coefs_h1) median(bias_coefs_h2) median(bias_coefs_h3) median(bias_coefs_h4)]
[median(ar1_coefficients),median(ar2_coefficients),median(ar3_coefficients),median(ar4_coefficients)]
[median(mz_coefbeta_g0_h1) median(mz_coefbeta_g0_h2) median(mz_coefbeta_g0_h3) median(mz_coefbeta_g0_h4)] 
[median(cg_h1_coefs) median(cg_h2_coefs) median(cg_h3_coefs)]