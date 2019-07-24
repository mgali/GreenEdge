% Preprocess CTD data and put in single mat file
clear, clc

%% Loop through casts

genpath = '~/Desktop/GreenEdge/CTD_rosette';
headrow = 36;
startrow = 38;

for ic = 1:203;
    
    fpath = sprintf('%s/1601_%03i.int',genpath,ic);
    
    if exist(fpath,'file')
        fid = fopen(fpath);
        formatDUMMY = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
        DUMMY = textscan(fid,formatDUMMY);
        fclose(fid);
        
        TMP = [];
        for j = 1:size(DUMMY,2)
            
            DUMMY{1,j}{36,1};
            tmp = DUMMY{1,j};
            head{j} = tmp{headrow};
            append = [];
            for k = startrow:size(tmp,1)
                append = [append; str2double(tmp{k})];
            end
            TMP = [TMP append];
        end
        
        ctd.DATA = TMP;
        ctd.head = head;
        CTD{ic} = ctd;
        
    else
        
        CTD{ic} = NaN;
        
    end
end

save(sprintf('%s/GreenEdge_CTD_2016_MGT.mat',genpath),'CTD')