% OPEN ANY ASCII (TEXT) FILE WITH FLEXIBLE SPECIFICATIONS
% Adapted from Anoop Mahajan's code and its application for my
% meta-analysis, Binning_conditioning.m, on 20170331

% % Comme once debugged
% filename = '~/Desktop/GreenEdge/Irradiance/COPS_profiles/GE2015.Ice.Camp/2015163_ICEPRO_150612_1419_C_data_004.tsv.asc/PAR.txt';
% headerline = 1;
% skiplines = 0;
% col_delimiter = '\t';
% NA_to_NaN = 1;
% Empty_to_NaN = 1;

function [head,DATA] = f_read_txt_dummy(filename,headerline,skiplines,col_delimiter,NA_to_NaN,Empty_to_NaN)

fid = fopen(filename);
DUMMY = textscan(fid,'%s','delimiter','\n'); % each row into a single string
fclose(fid); % close the file

% Extract header if any
if headerline
    head = regexp(DUMMY{1,1}{headerline,1},col_delimiter,'split'); % Create headers
else
    head = 'nan';
end

% Extract data
i = 0; % initialize row index
for j = (skiplines+1):numel(DUMMY{1,1})
    
    % Splitting the arrays into strings, comma delimited
    segmDATA = regexp(DUMMY{1,1}{j,1},col_delimiter,'split');
    i = i+1 ;
    for k = 1:numel(head),
        DATA{i,k} = segmDATA{1,k};
        
        % Replace NA with NaN?
        if NA_to_NaN
            if strcmp(DATA{i,k},'NA'), DATA{i,k} = 'NaN'; end
        end
        if Empty_to_NaN
            if isempty(DATA{i,k}), DATA{i,k} = 'NaN'; end
        end
    end
end