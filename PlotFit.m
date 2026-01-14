function  PlotFit(settings,tseries,fits,variances,model,state,titlestr)
% plotting function for state space results

if isempty(variances)
    PlotVariances=0;
else
    PlotVariances=1;
end

   % check existence irregular term
    if ~isfield(state,'epsilon') 
        disp('calculate residual')
        for p=1:tseries.ntseries
        % define as residual
        state.epsilon(p,:)=tseries.Y(p,:)-fits.fit{p};
        end
    end
    
 


for p=1:tseries.ntseries

    if settings.plotconstrainperiod
     % get indexes of period with data, including gaps   
      IndexData=[tseries.indexavailable(1):tseries.indexavailable(end)];
    else
      IndexData=[1:tseries.ntimes];
    end
    clear h
    % make new figure
    fig=figure;subplot(2,1,1);
    
    % get colors
    ColorOrder=get(gca,'ColorOrder');
    nColors=length(ColorOrder);
    
    for i=1:model.dim+1 % extra dimension for fit
        Colors(i,:)=ColorOrder(mod(i-1,nColors)+1,:);
    end
    
    % get index for which state is available, for plotting variance
    %if ~settings.plotconstrainperiod
        IndexDataVar=find(~isnan(fits.fit{p}));
    %end
   

    hold on
    h(1)=plot(tseries.time(IndexData),tseries.Y(p,IndexData),'k-');
    % plot without gaps
    plot(tseries.time(~isnan(tseries.Y(p,:))),tseries.Y(p,(~isnan(tseries.Y(p,:)))),'k:')
    plot(tseries.time,tseries.Y(p,:),'k.','MarkerSize',5);
    h(2)=plot(tseries.time(IndexData),fits.fit{p}(IndexData),'Color',Colors(1,:));
    
    h(3)=plot(tseries.time(IndexData),fits.trend{p}(IndexData),'Color',Colors(2,:));
    legendstring={'observation','fit','trend (1 \sigma)'};
    hi=3;
    if settings.intervention
        for i=1:settings.numberinterventions
            hi=hi+1;
            h(hi)=plot(tseries.time(IndexData),fits.intervention{p}(i,(IndexData)),'Color',Colors(hi,:));
            legendstring=[legendstring,{'step intervention'}];
            colorintervention=hi;
        end
    end
    
    if settings.regression
        for i=1:settings.numberregressors
            hi=hi+1;
            h(hi)=plot(tseries.time(IndexData),fits.regression{p}(i,(IndexData)),'Color',Colors(hi,:));
            legendstring=[legendstring,{'explanatory variable'}];
        end
    end
    
    
    if ~isfield(tseries,'name')
        tseries.name='';
    end
    if tseries.ntseries == 1
        if iscell(tseries.name)
            NamePlot=tseries.name{1};
        else
            NamePlot=tseries.name;
        end
    else
        NamePlot=strcat({tseries.name{p}});
    end
    
    xlabel('time')
    ylabel('total fit and trend')
    TitleString=strcat({titlestr},NamePlot);
    title(TitleString,'FontWeight','normal')
    box on
    
    if PlotVariances
        % plot 1-sigma values
        
        transparanterrorbars(tseries.time(IndexDataVar),fits.trend{p}(IndexDataVar),sqrt(variances.level{p}(IndexDataVar)),Colors(2,:));
        if settings.intervention
        for i=1:settings.numberinterventions
             variancesvec=ones(size(tseries.time(IndexDataVar)))'*sqrt(variances.delta{p}).*(fits.intervention{p}(IndexDataVar)~=0);
           transparanterrorbars(tseries.time(IndexDataVar),fits.intervention{p}(IndexDataVar),variancesvec,Colors(colorintervention,:));
        end
        end
    end
    
    legend(h,legendstring,'Location','EastOutside');
   % xlim([tseries.time(1) tseries.time(end)])
    xlim([tseries.time(IndexData(1)) tseries.time(IndexData(end))])
    subplot(2,1,2);hold on
    
 
    
    % plot
    clear  h
    hi=1;
    h(hi)=plot(tseries.time(IndexData),tseries.Y(p,(IndexData)),'k');
        % plot without gaps
    plot(tseries.time(~isnan(tseries.Y(p,:))),tseries.Y(p,(~isnan(tseries.Y(p,:)))),'k:')
    hi=hi+1;
    h(hi)=plot(tseries.time(IndexData),state.epsilon(p,(IndexData)),'Color',Colors(hi-1,:));
    legendstr={'observation','residual'};
    if settings.AR
        hi=hi+1;
        h(hi)=plot(tseries.time(IndexData),fits.AR{p}(IndexData),'Color',Colors(hi-1,:));
        legendstr{hi}='AR';
    end
    
    if settings.cycle
        for icycle=1:settings.numbercycles
            hi=hi+1;
            h(hi)=plot(tseries.time(IndexData),fits.cycle{p}(icycle,(IndexData)),'Color',Colors(hi-1,:));
        legendstr{hi}=strcat('cycle [',num2str(settings.periodscycle(icycle)),']');
        end
    end
            
    
    if PlotVariances
        hi=1;
        % plot 1-sigma values irregular
        if tseries.ntseries==1
            StDevIrr=sqrt(squeeze(state.varepsilon(:,:)))';
        else
        StDevIrr=sqrt(squeeze(state.varepsilon(p,p,:)))';
        end
        hi=hi+1;
        transparanterrorbars(tseries.time(IndexDataVar),state.epsilon(p,(IndexDataVar)),StDevIrr(IndexDataVar),Colors(hi-1,:));
        
        if settings.AR
            hi=hi+1;
            StDevAR=sqrt(variances.AR{p});
            transparanterrorbars(tseries.time(IndexDataVar),fits.AR{p}(IndexDataVar),StDevAR(IndexDataVar),Colors(hi-1,:));
        end
        
        
    if settings.cycle
        for icycle=1:settings.numbercycles
            hi=hi+1;
            StDevCycle=sqrt(variances.cycle{p}(icycle,:));
            h(hi)=transparanterrorbars(tseries.time(IndexDataVar),fits.cycle{p}(icycle,(IndexDataVar)),StDevCycle(IndexDataVar),Colors(hi-1,:));
        end
    end
        %   plot(tseries.time,state.epsilon(p,:)+StDevIrr,':c',tseries.time,state.epsilon(p,:)-StDevIrr,':','Color',Colors(2,:))
    end
    
    box on
    xlabel('time')
    ylabel('stationary components and residual')
    legend(h,legendstr,'Location','EastOutside')
 %   xlim([tseries.time(1) tseries.time(end)])
    xlim([tseries.time(IndexData(1)) tseries.time(IndexData(end))])
    % save figure
    if settings.saveplots
        
        
      if isfield(settings,'maindir')

          filename=char(strcat(settings.maindir,'/',settings.savedir,'/',{titlestr},char(NamePlot)));
      else
          filename=char(strcat(settings.savedir,'/',{titlestr},char(NamePlot)));
      end

        
        savefig(fig,filename);
        if settings.savepng
            print(fig,strcat(filename,'.png'),'-dpng')
        end
    end
end

end

