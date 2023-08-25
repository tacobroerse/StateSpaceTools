function [psd]=SpectralPlot(settings,tseries,estsmooth,fitsmooth);
% plot power spectral density

if settings.continuoustime || tseries.datagaps ~= 0

    


    % do Lomb-Scargle
    % loop on time series
    for p=1:tseries.ntseries

        % check for data gaps at beginning or end
        if isfield(tseries,'indexavailable')
            if tseries.ntseries==1
                IndexData=tseries.indexavailable;
            else
                IndexData=tseries.indexavailable{p};
            end
        else
            IndexData=[1:tseries.ntimes];
        end

    [psd{p}.Y,psd{p}.f,Palpha]=fastlomb(tseries.Y(p,(IndexData)),tseries.time(IndexData));
    [psd{p}.fitsmooth,dummy]=fastlomb(fitsmooth.fit{p}(IndexData),tseries.time(IndexData));
    [psd{p}.trendsmooth,dummy]=fastlomb(fitsmooth.trend{p}(IndexData),tseries.time(IndexData));
    [psd{p}.epsilonsmooth,dummy]=fastlomb(estsmooth.epsilon(p,(IndexData)),tseries.time(IndexData));
    if settings.cycle
        for i=1:settings.numbercycles
            [psd{p}.cycle{i},dummy]=fastlomb(fitsmooth.cycle{p}(i,(IndexData)),tseries.time(IndexData));
        end

    end

    if settings.AR
        [psd{p}.AR,dummy]=fastlomb(fitsmooth.AR{p}(IndexData),tseries.time(IndexData));
    end
    end


else
    for p=1:tseries.ntseries
    % do Fourier analysis
    [psd{p}.f,psd{p}.Y]=SpectralAnalysis(tseries.Y(p,:),tseries.time,tseries.frequency);
    [dummy,psd{p}.fitsmooth]=SpectralAnalysis(fitsmooth.fit{p},tseries.time,tseries.frequency);
    [dummy,psd{p}.trendsmooth]=SpectralAnalysis(fitsmooth.trend{p},tseries.time,tseries.frequency);
    [dummy,psd{p}.epsilonsmooth]=SpectralAnalysis(estsmooth.epsilon(p,:),tseries.time,tseries.frequency);
    if settings.cycle
        for i=1:settings.numbercycles
            [dummy,psd{p}.cycle{i}]=SpectralAnalysis(fitsmooth.cycle{i},tseries.time,tseries.frequency);
        end
    end
    end
end

 for p=1:tseries.ntseries
% Plot single-sided amplitude spectrum.
figure
legendstring={'observation','fit','trend','residual'};
semilogx(psd{p}.f,psd{p}.Y,'k') ; hold on
semilogx(psd{p}.f,psd{p}.fitsmooth)
semilogx(psd{p}.f,psd{p}.trendsmooth)
semilogx(psd{p}.f,psd{p}.epsilonsmooth)
if settings.cycle
    for i=1:settings.numbercycles
        legendstring=[legendstring,{strcat('cycle',num2str(i))}];
        semilogx(psd{p}.f,psd{p}.cycle{i});
    end
end

if settings.AR
legendstring=[legendstring,{'AR'}];
        semilogx(psd{p}.f,psd{p}.AR);
end

if tseries.ntseries==1
title(strcat({'Single-Sided Amplitude Spectrum of y(t) '},{tseries.name}),'FontWeight','normal')
else
   title(strcat({'Single-Sided Amplitude Spectrum of y(t) '},{tseries.name{p}}),'FontWeight','normal') 
end

xlabel('Frequency (1/yr)')
ylabel('|Y(f)|')
xlim([min(1/(tseries.ntimes*tseries.period)) max(psd{p}.f)])
legend(legendstring);
if settings.saveplots
    if isfield(tseries,'name')
        if tseries.ntseries==1
        filename=strcat(settings.savedir,'/','power_spectrum_',tseries.name);
        else
            filename=strcat(settings.savedir,'/','power_spectrum_',tseries.name{p});
        end
    else
        filename=strcat(settings.savedir,'/','power_spectrum');
    end
    savefig(filename);
end

 end
end

