% Test DMS and DMSPt algos with in situ GE2016 data
clear, clc

load ALLSURF
load samples.GE2016only.castsOnly.mat
GE2016 = samplesGE2016castsOnly;

% Choose dataset
% Surface data only
% load ALLSURF.20190213.48cols.mat
dmspt = ALLSURF(:,strcmp(headALLSURF,'dmspt'));
dms = ALLSURF(:,strcmp(headALLSURF,'dms'));
chl_in = ALLSURF(:,strcmp(headALLSURF,'Tchla'));

% % Data integrated over relevant layer
% % load ALLSURF.20190214.56cols.mat
% layer = 'z10';
% idms = ALLSURF(:,strcmp(headALLSURF,sprintf('idms_%s',layer)));
% idmspt = ALLSURF(:,strcmp(headALLSURF,sprintf('idmspt_%s',layer)));
% ichl = ALLSURF(:,strcmp(headALLSURF,sprintf('iTchla_%s',layer)));
% % zlayer = ALLSURF(:,strcmp(headALLSURF,'mld03'));
% zlayer = 10*ones(size(idms));
% dms = idms./zlayer;
% dmspt = idmspt./zlayer;
% chl_in = ichl./zlayer;

%% Calculate DMS-SAT DMS and DMSPt
zeutype = 'morel2007';
zeu_in = nan(size(chl_in)); 
if strcmp(zeutype,'morel2007')
    lchl = real(log10(chl_in));
    lzeu = 1.524 - 0.436*lchl - 0.0145*(lchl.^2) + 0.0186*(lchl.^3);
    zeu_in = 10.^lzeu;
end
mld_in = ALLSURF(:,strcmp(headALLSURF,'mld03'));
sst_in = ALLSURF(:,strcmp(headALLSURF,'sst'));
pic_in = nan(size(chl_in));
[dmspt_sat,flags_dsat,flag_key] = dmspt_algorithm_RSE2015(chl_in,zeu_in,mld_in,sst_in,pic_in);

paramDMS = [-1.300 0.700 0.0200];
par_in = ALLSURF(:,strcmp(headALLSURF,'PAR_at_3m_Einm-2d-1'));
par_in = 1.066*par_in/(exp(-0.06*3));
dms_sat = dms_sat(dmspt_sat,par_in,paramDMS);

exclude = (isnan(dms) & isnan(dmspt));
% exclude = exclude | (dms>10);
% exclude = exclude | (dmspt>100);
dms(exclude) = [];
dmspt(exclude) = [];
dms_sat(exclude) = [];
dmspt_sat(exclude) = [];
chl_in(exclude) = [];

ll = ''; % linear
skill_p.lin = f_skill_stats(dmspt, dmspt_sat, ll, '');
skill_d.lin = f_skill_stats(dms, dms_sat, ll, '');
skill_p_c.lin = f_skill_stats(dmspt, chl_in, ll, '');
skill_d_c.lin = f_skill_stats(dms, chl_in, ll, '');
[rlin.d, plin.d] = corr(dms, dms_sat,'rows','pairwise');
[rlin.p, plin.p] = corr(dmspt, dmspt_sat,'rows','pairwise');
[rs.d, ps.d] = corr(dms, dms_sat,'rows','pairwise','type','Spearman');
[rs.p, ps.p] = corr(dmspt, dmspt_sat,'rows','pairwise','type','Spearman');
ll = 'log10'; % log10, log (for ln)
skill_p.log = f_skill_stats(dmspt, dmspt_sat, ll, '');
skill_d.log = f_skill_stats(dms, dms_sat, ll, '');
skill_p_c.log = f_skill_stats(dmspt, chl_in, ll, '');
skill_d_c.log = f_skill_stats(dms, chl_in, ll, '');
[rlog.d, plog.d] = corr(real(log10(dms)), real(log10(dms_sat)),'rows','pairwise');
[rlog.p, plog.p] = corr(real(log10(dmspt)), real(log10(dmspt_sat)),'rows','pairwise');

%% Scatterplot DMSPt in situ vs SAT
figure(201), clf
loglog(dmspt, dmspt_sat,'.','markersize',40), hold on
l = legend(sprintf('N = %i',skill_p.log.N));
set(l,'location','northwest','fontsize',16,'box','on')
dm = nanmin([dmspt; dmspt_sat]);
dM = nanmax([dmspt; dmspt_sat]);
axis([dm dM dm dM])
plot([dm dM],[dm dM],'-k')
plot([dm dM],0.5*[dm dM],'--k')
plot([dm dM],2*[dm dM],'--k')
axis square
set(gca,'fontsize',12)
xlabel('DMSPt in situ (nM)','fontsize',16)
ylabel('DMSPt-SAT with in situ data input (nM)','fontsize',16)
title(sprintf('r_{log} = %0.2f, RMSE_{log} = %0.2f, bias_{lin} = %0.0f%%, MAPE_{lin} = %0.0f%%',...
    skill_p.log.r,skill_p.log.rms, 100*skill_p.lin.relbias, skill_p.lin.mape),'fontsize',16)

figure(202), clf
loglog(dms, dms_sat,'.','markersize',40), hold on
l = legend(sprintf('N = %i',skill_d.log.N));
set(l,'location','northwest','fontsize',16,'box','on')
dm = nanmin([dms; dms_sat]);
dM = nanmax([dms; dms_sat]);
axis([dm dM dm dM])
plot([dm dM],[dm dM],'-k')
plot([dm dM],0.5*[dm dM],'--k')
plot([dm dM],2*[dm dM],'--k')
axis square
set(gca,'fontsize',12)
xlabel('DMS in situ (nM)','fontsize',16)
ylabel('DMS-SAT with in situ data input (nM)','fontsize',16)
title(sprintf('r_{log} = %0.2f, RMSE_{log} = %0.2f, bias_{lin} = %0.0f%%, MAPE_{lin} = %0.0f%%',...
    skill_d.log.r,skill_d.log.rms, 100*skill_d.lin.relbias, skill_d.lin.mape),'fontsize',16)
