% PLOT TEST OF run_get_ed0.m

clear
clc

%% Load data

mainpath = '~/Desktop/GreenEdge/Irradiance/Ed_RT_MODISA/2016/280';
above = 0;
interpol = 0;
ed.plus = dlmread(sprintf('%s/interpol_%01i/Ed0plus_2016_280_69.325_-60.967.txt',mainpath,interpol));
ed.moins = dlmread(sprintf('%s/interpol_%01i/Ed0moins_2016_280_69.325_-60.967.txt',mainpath,interpol));

%% Plot

wl = 290:5:700;
hours = 0:3:21;
color.moins = [1 .3 .3];
color.plus = [0 .8 .8];

figure(99), clf
for ih = 1:size(ed.plus,2)
    
    subplot(3,3,ih)
    plot(wl,ed.moins(:,ih),'-','linewidth',2,'color',color.moins), hold on
    plot(wl,ed.plus(:,ih),'-','linewidth',2,'color',color.plus)
    axis([280 710 0 1]);
    set(gca,'tickdir','out')
    text(300,1.05,sprintf('Y=2016, DOY=280, %02i:00:00',hours(ih)),'fontsize',12)
    grid on
    
    if ih == size(ed.plus,2)
        text(900,0.6,'--- Ed0moins','color',color.moins,'fontsize',12)
        text(900,0.4,'--- Ed0plus','color',color.plus,'fontsize',12)
    end
    
end