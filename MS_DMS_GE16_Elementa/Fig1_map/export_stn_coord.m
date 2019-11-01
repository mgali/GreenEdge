% Export station coordinates list
clear, clc
load ~/Desktop/GreenEdge/GCMS/samples.GE2016only.castsOnly.mat
toexport = samplesGE2016castsOnly(:,[1 2 6 7 8]);
toexport(isnan(toexport(:,end)),:) = [];
csvwrite('GE_Amundsen_Station_Coordinates.csv',toexport)
toexport(toexport(:,1)<400,:) = [];
csvwrite('GE_Amundsen_Station_Coordinates_leg1b.csv',toexport)