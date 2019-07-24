% LOAD ALL MERGED *.mat FILES AND EXPORT AS *.csv
% NOTE: using array2table to format data and easily export matrices along
% with their headers stored in cell arrays. Otherwise, headers cannot be
% easily written using dlmwrite or similar.

clear
clc

%% List of stations and surface data
load samples.GE2016only.castsOnly.mat

% ALLSURF: master file for merged surface data
load ALLSURF

% Add cast number
cast = (1:203)';

% Matrix for export
M = [cast samplesGE2016castsOnly ALLSURF];

% Convert all variable names to suitable format
headALLSURF{2} = 'isolume_m_at_01';
headALLSURF{3} = 'isolume_m_at_0415';
headALLSURF{5} = 'MixedLayerDepth_m_01_criterion';
headALLSURF{8} = 'AW_ArW_clustering_coefficient';
headALLSURF{10} = 'PAR_at_3m_E_m2_d';

% Header for export
hM = ['cast'; headerGE2016; headALLSURF];

% Create table
tableforexport = array2table(M,'VariableNames',hM);

% Write
writetable(tableforexport,'GE2016.casts.ALLSURF.csv')

%% Profile data

% load surf_profile_ALL_4plot.20190123.mat
load surf_profile_ALL_4plot

% Delete rows labelled 'nouse'
rm = zeros(size(headDATA_CTD));
for j = 1:length(headDATA_CTD)
    if length(headDATA_CTD{j})==5 && sum(headDATA_CTD{j}=='nouse')==5, rm(j) = 1; end
end
DATA_CTD(:,rm==1) = [];
headDATA_CTD(rm==1) = [];

% Create table
tableforexport2 = array2table(DATA_CTD,'VariableNames',headDATA_CTD);

% Write
writetable(tableforexport2,'GE2016.profiles.DATA_CTD.DMS.csv')

