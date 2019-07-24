% PLOT DMS PROFILES FROM GREEN EDGE
clear
clc
data = xlsread('~/Desktop/GreenEdge/GCMS/gcms_greenedge_Marti.xlsx','Niskin_underway','A2:AC135');

% Sort
[data, sortindexs] = sortrows(data,[1 2 3]);

prescreen = data(:,15);
data(prescreen


