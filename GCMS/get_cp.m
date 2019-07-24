% Calculate beam attenuation (cp) from transmittance profiles.
% Make and save plots

function cp_out = get_cp(t_in, z_in, kstn)
% t_in = Trans; z_in = depth;

tmp = [z_in t_in];
exclude = sum(abs(tmp),2) == Inf | sum(isnan(tmp),2) > 0;
z_in(exclude) = [];
t_in(exclude) = [];

kstn = kstn*10; % deal with stn 604.5
pl = 0.25; % pathlength, [m]
zref = [150 200]; % reference depth horizon for baseline
c_w = 0.358; % cp_600 of pure water, [m-1]. After Loisel&Morel1998, Bishop1986.
t_inmax = nanmax(t_in);
t_in = 100*t_in/t_inmax;
c_raw = log(100./t_in)/pl; % raw c, [m-1]
zmax = nanmax(z_in);
if zmax >= zref(1)
    c_b = nanmean(c_raw(z_in>=zref(1)&z_in<=zref(2))); % baseline c: fixed depth horizon
    c_q = nanmean([nanmin(c_raw) quantile(c_raw,0.01)]); % baseline c: lowest cp in profile
else
    c_b = 0;
    c_q = 0;
end
cp_out.b = c_raw - c_b; % corrected for baseline c
cp_out.q = c_raw - c_q; % corrected for lowest c

% Smooth the minimum-corrected cp
cp_out.s = smooth(cp_out.q,5);
cp_out.ss = smooth(cp_out.s,5);

% Calculate 0.1 ad 0.9 quantiles of raw-smooth residuals
resid = cp_out.q./cp_out.ss;

% Identify positive spikes
qhigh = quantile(resid, .80);
ispikes = resid > qhigh;

% Despike
cp_out.d = cp_out.s;
cp_out.d(ispikes) = nan;
cp_out.d(ispikes) = interp1(z_in(~ispikes), cp_out.s(~ispikes), z_in(ispikes));

% Find depth of maximum cp
cp4zmax = cp_out.d;
cp4zmax(z_in<3 | z_in>50) = nan;
cp_out.zdbm = mode(z_in(cp4zmax==nanmax(cp4zmax)));

% Ensure output has same length as input
cp_out.b = [cp_out.b; nan(sum(exclude),1)];
cp_out.q = [cp_out.q; nan(sum(exclude),1)];
cp_out.s = [cp_out.s; nan(sum(exclude),1)];
cp_out.ss = [cp_out.ss; nan(sum(exclude),1)];
cp_out.d = [cp_out.d; nan(sum(exclude),1)];
 

iall = (1:length(exclude))';
reorder = [iall(~exclude); iall(exclude)];

cp_out.b = cp_out.b(reorder);
cp_out.q = cp_out.q(reorder);
cp_out.s = cp_out.s(reorder);
cp_out.ss = cp_out.ss(reorder);
cp_out.d = cp_out.d(reorder);

% % Plot and save
% figure(kstn), clf
% set(gcf,'units','centimeters','position',[2 2 40 16])
% 
% subplot(131)
% plot(t_in,-z_in,'-k')
% title(sprintf('Transmittance \n\n t-inmax = %0.2f',t_inmax))
% xlabel('T (%)')
% 
% subplot(132)
% plot(cp_out.b,-z_in,'-','color',[.6 .6 .6]), hold on
% plot(cp_out.q,-z_in,'-b')
% plot(cp_out.d,-z_in,'-m','linewidth',3)
% plot(cp_out.s,-z_in,'-c','linewidth',1)
% plot([0 0],-[nanmin(z_in) nanmax(z_in)],'-k')
% plot(cp_out.d(z_in==cp_out.zdbm), -cp_out.zdbm, 'ok', 'markersize', 10, 'linewidth', 3)
% title(sprintf('Corrected cp600 \n\n c_{bl} = %0.3f, c_{q} = %0.3f, ',c_b,c_q))
% xlim([-0.01 1.05*nanmax([cp_out.b; cp_out.q])])
% legend('150-200 m subtracted','min subtracted','min sub+smooth+despiked','min sub+smoothed','location','south')
% xlabel('cp600 (m^{-1})')
% 
% subplot(133)
% hist(resid, 100), hold on
% line(ones(1,2)*qhigh,[0 200],'color','r','linewidth',1)
% set(gca,'xscale','log','xlim',[0.3 10])
% 
% dirout = '~/Desktop/GreenEdge/GCMS/plots_cp_profile_v201901';
% print(kstn,sprintf('%s/Trans_cp_%i.png',dirout,kstn),'-dpng','-r300')
% close(kstn)
