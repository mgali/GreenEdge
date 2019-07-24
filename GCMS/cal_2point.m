% CALCULATE DMS 62 FROM 2-POINT REGRESSION BETWENN 65-68 PA AND nM DMS
function dms62 = cal_2point(DATA,correct_frozen,corr_pa)

dms62 = nan(size(DATA,1),1);

for j = 1:size(DATA,1)
    
    injvol = DATA(j,14); % injection volume in mL
    pa62 = DATA(j,16); % peak areas
    pa65 = DATA(j,17);
    pa68 = DATA(j,18);
    ec65 = DATA(j,25); % expected concentrations in nmol/L
    ec68 = DATA(j,26);
    
    if pa62 >= pa65 && pa62 <= pa68
        X = [1 pa65; 1 pa68];
        y = [ec65; ec68];
        fit = X\y;
        correctnM = (fit(1) + fit(2)*corr_pa)*injvol/20;
        dms62(j) = (fit(1) + fit(2)*pa62)*injvol/20 - correctnM;
    end
end

dms62 = dms62.*correct_frozen;