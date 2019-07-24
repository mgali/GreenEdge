% READ XLS WITH TAXONOMY, MICROSCOPY COUNTS 
% MARTI GALI TAPIAS 20180705

function [TAXO,headTAXO] = f_readxls_taxo(filepath,filename)

headTAXO = {'cast' 'depth' 'diat_cen' 'diat_pen' 'dino_athec' 'dino_thec'...
    'chlor' 'chrys' 'dictyo' 'crypt' 'eugl' 'pras' 'prym' 'Phaeo' 'flag'...
    'raph' 'cyan' 'hetero' 'choano' 'cilli'}'; % 'unID_2_5' 'unID_5_10' 'unID_10_20'
col = {'L' 'N' 'T' 'BQ' 'DJ' 'ED'...
    'FB' 'FC' 'FK' 'FN' 'FT' 'FX' 'GG' 'GL' 'GP'...
    'GV' 'GX' 'GZ' 'HC' 'HK'}'; %  'IL' 'IM' 'IN'
TAXO = nan(length(2:37),length(headTAXO));
if length(headTAXO)==length(col);
    for j = 1:length(headTAXO)
        % sprintf('column %i, xls range %s2:%s37',j,col{j},col{j})
        tmp = xlsread(sprintf('%s/%s',filepath,filename),...
            sprintf('%s2:%s37',col{j},col{j}));
        % sprintf('matrix dims n=%i, m=%i',size(tmp,1),size(tmp,2))
        TAXO(:,j) = tmp;
    end
else
    error('Check headers and matching column indices')
end

% diatoms = S;
% centric = T;
% pennate = BQ;
% dino_athec = DJ;
% dino_thec = ED;
% chloro = FB;
% chryso = FC;
% dictyo = FK;
% crypto = FN;
% eugle = FT;
% prasi = FX;
% prym = GG;
% flag = GP;
% raph = GV;
% cyano = GX;
% hetero = GZ;
% choano = HC;
% cilli = HK;
% unID_2_5 = IL;
% unID_5_10 = IM;
% unID_10_20 = IN;
