% Calculate MLD from surface reference depth = 2m

function [mld,N2max,zN2max] = get_mld(depth,sigt,N2,kstn)

kstn = kstn*10;

% Smooth
orig.sigt = sigt;
orig.N2 = N2;
sigt = smooth(sigt,5);
N2 = smooth(N2,5);

sref = sigt(depth<2);
if isempty(sref),
    sref = nanmin(depth(depth<4));
end
i03 = 0;
dsigt = 0;
while dsigt < 0.03
    i03 = i03+1;
    dsigt = abs(sigt(i03) - sref);
    mld.d03 = depth(i03);
end
i125 = 0;
dsigt = 0;
while dsigt < 0.125
    i125 = i125+1;
    dsigt = abs(sigt(i125) - sref);
    mld.d125 = depth(i125);
end

N2d03 = N2(depth<=mld.d125 & depth>=mld.d03);
zd03 = depth(depth<=mld.d125 & depth>=mld.d03);
N2max.d03 = nanmax(N2d03);
if ~isnan(N2max.d03), zN2max.d03 = zd03(N2max.d03==N2d03); else zN2max.d03 = nan; end
N2d125 = N2(depth<=100 & depth>=mld.d125);
zd125 = depth(depth<=100 & depth>=mld.d125);
N2max.d125 = nanmax(N2d125);
if ~isnan(N2max.d125), zN2max.d125 = zd125(N2max.d125==N2d125); else zN2max.d125 = nan; end

% % Plot and save
% figure(kstn), clf
% subplot(121)
% plot(orig.sigt,-depth,'-','linewidth',1,'color',[.6 .6 .6]), hold on
% plot(sigt,-depth,'-k','linewidth',2)
% plot(sigt(i03),-depth(i03),'.c','markersize',30)
% plot(sigt(i125),-depth(i125),'.b','markersize',30)
% ylim([-100 0])
% xlabel('sigma-t (kg m^{-3})')
% ylabel('Depth')
% legend('0.2 m bins','1m running mean','MLD03','MLD125','location','south')
% 
% subplot(122)
% plot(orig.N2,-depth,'-','linewidth',1,'color',[.6 .6 .6]), hold on
% plot(N2,-depth,'-k','linewidth',2), hold on
% plot(N2max.d03,-zN2max.d03,'.y','markersize',30)
% plot(N2max.d125,-zN2max.d125,'.r','markersize',30)
% ylim([-100 0])
% xlabel('d(sigma-t)/dz (kg m^{-4})')
% ylabel('Depth')
% legend('0.2 m bins','1m running mean','N2-MLD03','N2-MLD125','location','south')
% 
% print(kstn,sprintf('Sigt_MLD_N2_%i.png',kstn),'-dpng','-r300')
% close(kstn)





