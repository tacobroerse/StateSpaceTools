function  PlotMeanSlope(tseries,settings,fitsmooth,meanslope,varsmooth,model)
% plot mean slope estimated from state space results

if settings.estimatemeanslope
    if settings.slope
        
        fig=figure;hold on
        box on
        
               % get colors
    ColorOrder=get(gca,'ColorOrder');
    nColors=length(ColorOrder);
    j=0;
    for p=1:tseries.ntseries
    for i=1:model.dim+1 % extra dimension for fit
        j=j+1;
        Colors(j,:)=ColorOrder(mod(i-1,nColors)+1,:);
    end
    end  
    
    jPlot=1;
        for p=1:tseries.ntseries
        % get averaged slope and error variance
        if tseries.ntseries==1
            name=tseries.name;
        else
        name=tseries.name{p};
        end
        % 
   
        
%         if settings.lsq
%             % fit a least squares slope to the time variable slope
%             tseriestemp=tseries;
%             tseriestemp.y=fitsmooth.trend';
%             settingstemp=settings;
%             settingstemp.cycle=0;
%             [estlsqtrend]=LeastSquaresTrend(tseriestemp,settingstemp,p);
%         end
     
        legendstr{jPlot}=['slope ',name];
        transparanterrorbars(tseries.time,fitsmooth.slope{p},sqrt(varsmooth.slope{p}),ColorOrder(jPlot,:));
        h(jPlot)=plot(tseries.time,fitsmooth.slope{p},'Color',ColorOrder(jPlot,:));jPlot=jPlot+1;
         
        legendstr{jPlot}=['mean slope ',name];
         transparanterrorbars(tseries.time,ones(1,tseries.ntimes)*meanslope.slope{p},ones(1,tseries.ntimes)*sqrt(meanslope.slopevar{p}),ColorOrder(jPlot,:));
        h(jPlot)=plot(tseries.time,ones(1,tseries.ntimes)*meanslope.slope{p},'Color',ColorOrder(jPlot,:));jPlot=jPlot+1;
        
        
%         if settings.lsq
%             plot(tseries.time,ones(1,tseries.ntimes)*estlsqtrend.slope)
%         end
%         
%         if settings.lsq
%             
%             plot(tseries.time,ones(1,tseries.ntimes)*estlsqtrend.slope,'-r')
%             
%             plot(tseries.time,ones(1,tseries.ntimes)*(estlsqtrend.slope+sqrt(estlsqtrend.cov(2,2))),':r')
%             plot(tseries.time,ones(1,tseries.ntimes)*(estlsqtrend.slope-sqrt(estlsqtrend.cov(2,2))),':r')
%            
%         end
       
        
        
        
        
%         if settings.lsq
%             legend('slope','mean slope','least squares slope')
%             if settings.lsq
%                 legend('slope','mean slope','lsq trend','least squares slope')
%             end
%         end
         
        end % end loop on time series
        
        xlabel('time')
        ylabel('slope')
        legend(h,legendstr)
        
        title('slope and mean slope and 1-\sigma errors ','FontWeight','normal')
        % save figure
        if settings.saveplots
            if isfield(tseries,'name')
                
                filename=strcat(settings.maindir,'/',settings.savedir,'/','mean_slope_',tseries.generalname);
            else
                filename=strcat(settings.maindir,'/',settings.savedir,'/','mean_slope');
            end
            savefig(filename);
            if settings.savepng
                print(fig,strcat(filename,'.png'),'-dpng')
            end
        end
    end
end

end

