function PlotCovariance(tseries,model,covsmooth,estsmooth,settings,typeplot)
% PlotCovariance plots covariance matrices for the estimated state
%
% HOW:     PlotCovariance(tseries,model,covsmooth,estsmooth,settings,typeplot)
%
% Input:
%       tseries structure which should contain at least:
%           ntimesobs                      number of epochs
%
%       model                          contains indices of components
%
%       covsmooth                       contains the correlation between
%           states at different epochs
%
%       estsmooth                      V constains variances and
%           covariances per epoch
%
%       typeplot                    choose 'autocorrelation',
%       'crosstimeseriescorrelation','crosscomponentcorrelation'
%
%
%
% Note:
%
% Taco Broerse, Utrecht University 2020
% d.b.t.broerse@uu.nl

% colors
% colors from scientific colormap, roma has large range of colors
%colormaps=load('roma.mat');
colormaps=load('cork.mat');
% roma=colormaps.roma;
% cmap=flipud(roma);
% colormaps=load('cork.mat');
cmap=colormaps.cork;




switch typeplot
    
    case 'autocorrelation'
        % auto correlation (i.e. correlation of a state with itself at a
        % different epoch)
        
        % set number of panels
        n=tseries.ntseries;
        m=1;
        if settings.slope
            m=m+1;
        end
        if settings.AR
            m=m+1;
        end
        
        % set correlation matrices
        
        for p=1:tseries.ntseries
            if tseries.ntseries==1
                tseriesname=tseries.name;
            else
                tseriesname=tseries.name{p};
            end
            % trend level
            % correlation from covariance
            q=1;
            CorMat{p,q}=squeeze(covsmooth.corr_coef(model.indexlevel{p},model.indexlevel{p},:,:));
            % correlation 1 on the diagonal
            for t=1:tseries.ntimesobs
                CorMat{p}(t,t)=1;
            end
            
            
            % trend slope
            if settings.slope
                q=q+1;
                % correlation from covariance
                CorMat{p,q}=squeeze(covsmooth.corr_coef(model.indexslope{p},model.indexslope{p},:,:));
                % correlation 1 on the diagonal
                for t=1:tseries.ntimesobs
                    CorMat{p,q}(t,t)=1;
                end
            end
            
            % AR
            if settings.AR
                q=q+1;
                % correlation from covariance
                CorMat{p,q}=squeeze(covsmooth.corr_coef(model.indexAR{p},model.indexAR{p},:,:));
                % correlation 1 on the diagonal
                for t=1:tseries.ntimesobs
                    CorMat{p,q}(t,t)=1;
                end
            end
            
            % set titles
            maintitle=['state autocorrelation ' tseries.generalname];
            q=1;
            titlestr{p,q}=[tseriesname newline 'trend level'];
            xlabelstr{p,q}='time';
            ylabelstr{p,q}='time';
            if settings.slope
                q=q+1;
                titlestr{p,q}=[tseriesname newline 'trend slope'];
                xlabelstr{p,q}='time';
                ylabelstr{p,q}='time';
            end
            if settings.AR
                q=q+1;
                titlestr{p,q}=[tseriesname newline 'AR'];
                xlabelstr{p,q}='time';
                ylabelstr{p,q}='time';
            end
            
            
            
            
            
            
        end
        
    case 'crosstimeseriescorrelation'
        
        % correlation between different components, different time series
        
        if tseries.ntseries==1
            error('only possible for multivariate models')
        end
        % set number of panels
        
        if tseries.ntseries> 2
            error('not implemented for more than two time series')
        end
        m=1;
        if settings.slope
           m=m+1; 
        end
        if settings.AR
            m=m+1;
        end
        
        n=m;
        
        for p=1:n
            if tseries.ntseries==1
                tseries.tseriesname{p}=tseries.name;
            else
                tseries.tseriesname{p}=tseries.name{p};
            end
            
            for q=1:m
                if p==1
                    index1=model.indexlevel{1};
                elseif p==2
                    index1=model.indexslope{1};
                elseif p==3
                    index1=model.indexAR{1};
                end
                if q==1
                    index2=model.indexlevel{2};
                elseif q==2
                    index2=model.indexslope{2};
                elseif q==3
                    index2=model.indexAR{2};
                end
                
                CorMat{p,q}=squeeze(covsmooth.corr_coef(index1,index2,:,:));
                for t=1:tseries.ntimesobs
                    CorMat{p,q}(t,t)=estsmooth.V(index1,index2,t)...
                        /sqrt(estsmooth.V(index1,index1,t)*...
                        estsmooth.V(index2,index2,t));
                end
            end
        end
        
      
       
        
        % set titles
        maintitle=['state component crosscorrelation '  tseries.generalname];
        for p=1:n
            for q=1:m
                if p==1
                    type1='trend ';
                    ylabelstr{p,q}=['time ' tseries.tseriesname{1}];
                elseif p==2
                    type1='slope ';
                    ylabelstr{p,q}=['time ' tseries.tseriesname{1}];
                elseif p==3
                    type1='AR ';
                    ylabelstr{p,q}=['time ' tseries.tseriesname{1}];
                end
                if q==1
                    type2='trend ';
                    xlabelstr{p,q}=['time ' tseries.tseriesname{2}];
                elseif q==2
                    type2='slope ';
                    xlabelstr{p,q}=['time ' tseries.tseriesname{2}];
                elseif q==3
                    type2='AR ';
                    xlabelstr{p,q}=['time ' tseries.tseriesname{2}];
                end
                titlestr{p,q}=[type1 tseries.tseriesname{1} ' vs. ' type2 tseries.tseriesname{2}];
            end
        end
        
      
        
        % end
        
    case 'crosscomponentcorrelation'
        % correlation between different components
        
        
        % set number of panels
        n=tseries.ntseries;
        m=1;
        
        if settings.AR
            m=m+1;
        end
        
        for p=1:tseries.ntseries
            if tseries.ntseries==1
                tseries.tseriesname{p}=tseries.name;
            else
                tseries.tseriesname{p}=tseries.name{p};
            end
            
            % correlation between trend level and trend slope
            q=1;
            CorMat{p,q}=squeeze(covsmooth.corr_coef(model.indexlevel{p},model.indexslope{p},:,:));
            for t=1:tseries.ntimesobs
                CorMat{p,q}(t,t)=estsmooth.V(model.indexlevel{p},model.indexslope{p},t)...
                    /sqrt(estsmooth.V(model.indexlevel{p},model.indexlevel{p},t)*estsmooth.V(model.indexslope{p},model.indexslope{p},t));
            end
            % correlation between trend (level) and AR
            if settings.AR
                q=q+1;
                CorMat{p,q}=squeeze(covsmooth.corr_coef(model.indexlevel{p},model.indexAR{p},:,:));
                for t=1:tseries.ntimesobs
                    CorMat{p,q}(t,t)=estsmooth.V(model.indexlevel{p},model.indexAR{p},t)...
                        /sqrt(estsmooth.V(model.indexlevel{p},model.indexlevel{p},t)*estsmooth.V(model.indexAR{p},model.indexAR{p},t));
                end
            end
            
            
            % set titles
            maintitle=['state component crosscorrelation ' tseries.generalname];
            q=1;
            titlestr{p,q}=[tseries.tseriesname{p} newline 'trend level vs trend slope'];
                xlabelstr{p,q}='time trend slope';
                ylabelstr{p,q}='time trend level';
            
            if settings.AR
                q=q+1;
                titlestr{p,q}=[tseries.tseriesname{p} newline 'AR vs. trend'];
                xlabelstr{p,q}='time AR';
                ylabelstr{p,q}='time trend level';
            end
        end
        
    otherwise
        
        error(strcat('invalid option for typeplot :',typeplot))
        
end

% make figure
FigurePanelHeight=300;
fig=figure;
fig.Position=[100 100 FigurePanelHeight*m FigurePanelHeight*n];
colormap(cmap)
for p=1:n
    for q=1:m
        subplot(n,m,(p-1)*m+q);
        imagesc(tseries.time,tseries.time,CorMat{p,q});colorbar
        title(titlestr{p,q},'FontWeight','normal')
        axis square
        caxis([-1 1])
        axisi=gca;
        % place x labels on top
        axisi.XAxisLocation='bottom';
        axisi.YTickLabel=axisi.XTickLabel;
        axisi.YTick=axisi.XTick;
        if p==n
           xlabel(xlabelstr{p,q})
        end
        if q==1
            ylabel(ylabelstr{p,q})
        end
       % end
    end
    
    
end
sgtitle(maintitle,'FontWeight','normal')
% save figure
if settings.saveplots
    if isfield(tseries,'name')
        filename=strcat(settings.savedir,'/',typeplot,'_',tseries.generalname);
    else
        filename=strcat(settings.savedir,'/',typeplot);
    end
    savefig(filename);
    if settings.savepng
        print(fig,strcat(filename,'.png'),'-dpng')
    end
end

end

