% Export station coordinates list
clear, clc
load ~/Desktop/GreenEdge/GCMS/samples.GE2016only.castsOnly.mat
load ~/Desktop/GreenEdge/GCMS/ALLSURF
toexport = [samplesGE2016castsOnly(:,[1 2 8 6 7]) ALLSURF(:,[13 31 32 33])];
toexport(isnan(toexport(:,end)),:) = [];
csvwrite('GE_Amundsen_fdms_sicday.csv',toexport)
