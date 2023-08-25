function  PlotCycles(settings,tseries,fits,variances)
% plotting function for state space results, cycles only
% plot cycles and their amplitudes

if settings.cycle
    for p=1:tseries.ntseries
        fig=figure;clear h
        
if settings.plotconstrainperiod
     % get indexes of period with data, including gaps   
      IndexData=[tseries.indexavailable(1):tseries.indexavailable(end)];
    else
      IndexData=[1:tseries.ntimes];
    end

      % get index for which state is available, for plotting variance
    %if ~settings.plotconstrainperiod
        IndexDataVar=find(~isnan(fits.fit{p}));
    %end
        ax1=subplot(2,1,1);
        ax1Position=ax1.Position;
        hold on
        % get colors
        ColorOrder=get(gca,'ColorOrder');
        
        
        legendstring=[];
        
        for i=1:settings.numbercycles
            
            h(i*2-1)=plot(tseries.time,fits.cycle{p}(i,:),'-','Color',ColorOrder(i,:));
            
            h(i*2)=plot(tseries.time,fits.cycleampl{p}(i,:),':','Color',ColorOrder(i,:));
            legendstring=[legendstring,{strcat('cycle with period:',num2str(settings.periodscycle(i)))}];
            legendstring=[legendstring,{'amplitude'}];
        end
        if tseries.ntseries==1
            title(strcat({'cycle components '},tseries.name),'FontWeight','normal')
        else
        title(strcat({'cycle components '},{tseries.name{p}}),'FontWeight','normal')
        end
        
        ylabel('amplitude')
        
        
        % plot 1-sigma values
        for i=1:settings.numbercycles
            StDevCycle=sqrt(variances.cycle{p}(i,:));
            transparanterrorbars(tseries.time(IndexDataVar),fits.cycle{p}(i,IndexDataVar),StDevCycle(IndexDataVar),ColorOrder(i,:));
            
        end
      %  xlim([tseries.time(1) tseries.time(end)])
       xlim([tseries.time(IndexData(1)) tseries.time(IndexData(end))])
        legend(h,legendstring,'Location','Best','Box','off')
        ax1.Position=ax1Position;
        
        
        % plot phases
        ax2=subplot(2,1,2);
        
        hold on
        for i=1:settings.numbercycles
            plot(tseries.time,fits.cyclephase{p}(i,:),'.','Color',ColorOrder(i,:))
            
        end
        title('cycle components phase offset')
        ylim([0 2*pi])
        xlabel('time')
        ylabel('phase [rad]')
      %  xlim([tseries.time(1) tseries.time(end)])
         xlim([tseries.time(IndexData(1)) tseries.time(IndexData(end))])
        % save figure
        if settings.saveplots
            if isfield(tseries,'name')
                if tseries.ntseries==1
                    filename=strcat(settings.maindir,'/',settings.savedir,'/','cycles_',tseries.name);
                else
                filename=strcat(settings.maindir,'/',settings.savedir,'/','cycles_',tseries.name{p});
                end
            else
                filename=strcat(settings.maindir,'/',settings.savedir,'/','cycles');
            end
            savefig(filename);
            if settings.savepng
                print(fig,strcat(filename,'.png'),'-dpng')
            end
        end
    end
end



end

