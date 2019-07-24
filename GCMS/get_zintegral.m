% Calculate beam attenuation (cp) from transmittance profiles.
% Make and save plots

function varint_out = get_zintegral(z_in, var_in, zbounds, kk, kcast, varname, layername)
if length(z_in) == length(var_in)
    kk = kk*10;
    profile = [z_in var_in];
    % Remove nan-containing rows and data outside bounds
    exclude = z_in<zbounds(1) | z_in>zbounds(2) | sum(isnan(profile),2)~=0;
    profile(exclude,:) = [];
    profile = sortrows(profile,1);
    % Average together measurements repeated at same depth
    zdiff = diff(profile(:,1));
    if sum(zdiff==0)
        znew = [profile(zdiff~=0,1); profile(end,1)];
        profilenew = [];
        for k = 1:length(znew)
            tmp = profile(profile(:,1)==znew(k),:);
            profilenew = [profilenew; nanmean(tmp,1)];
        end
        profile = profilenew;
    end
    % Check n of measurements in layer being integrated (see *Alternative)
    % Minimum is 1 measurement every 11 m for layers <=40 m and 1
    % measurement every 16 m for layers > 40 m
    nmeas = sum(~isnan(profile(:,2)));
    if (zbounds(2)/nmeas <= 11 && zbounds(2) <= 40) || (zbounds(2)/nmeas <= 16 && zbounds(2) > 40)
        % Perform linear vertical inerpolation
        z_out = (1:100)';
        profile = [profile; [100 0]]; % put 0 nM at 100 m
        var_out = interp1(profile(:,1), profile(:,2), z_out);
        varint_out = nanmean(var_out((zbounds(1)+1):zbounds(2)),1)*(diff(zbounds));
        % Make profile plots
        if zbounds(2)>=60
            figure(kcast)
            plot(profile(:,2), -profile(:,1), '.r', 'markersize', 20), hold on
            plot(var_out, -z_out, '.-k', 'markersize', 6)
            legend('original','interpolated','fontsize',16,'location','southeast')
            plot([0 nanmax(var_out)],-zbounds(1)*ones(1,2),'--b','linewidth',2)
            plot([0 nanmax(var_out)],-zbounds(2)*ones(1,2),'--b','linewidth',2)
            title(sprintf('Vertical integral = %0.1f units m^{-2}',varint_out),...
                'fontsize',20)
            dirout = '~/Desktop/GreenEdge/GCMS/plots_zinterp';
            print(kcast,sprintf('%s/profile_interp_st%i_cast%i_%s_%s.png',dirout,kk,kcast,varname,layername),'-dpng','-r300')
            close(kcast)
        end
    else
        varint_out = nan;
    end
else
    varint_out = nan;
    warning('Non-matching vector lengths, no output produced')
end

% *Alternative: Minimum is 2 for layer ?20 m and 4 for layer >20 m
% if (nmeas >= 2 && zbounds(2)<=20) || (nmeas >= 4 && zbounds(2)>20)
