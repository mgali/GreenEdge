% Calculate beam attenuation (cp) from transmittance profiles.
% Make and save plots

function cp_out = get_cp(t_in, z_in, kstn)
% t_in = Trans; z_in = depth;

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

% Plot and save
figure(kstn), clf
subplot(121)
plot(t_in,-z_in,'-k')
title(sprintf('Corrected transmittance \n\n t-inmax = %0.2f',t_inmax))
xlabel('T (%)')

subplot(122)
plot(cp_out.b,-z_in,'-r'), hold on
plot(cp_out.q,-z_in,'-b')
plot(cp_out.s,-z_in,'-c','linewidth',2)
plot([0 0],-[nanmin(z_in) nanmax(z_in)],'-k')
title(sprintf('Corrected cp600 \n\n c_{bl} = %0.3f, c_{q} = %0.3f, ',c_b,c_q))
xlim([-0.01 1.05*nanmax([cp_out.b; cp_out.q])])
legend('150-200 m subtracted','min subtracted','min sub., smoothed','location','south')
xlabel('cp600 (m^{-1})')

print(kstn,sprintf('Trans_cp_%i.png',kstn),'-dpng','-r300')
close(kstn)
