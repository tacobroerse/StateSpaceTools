function [distvar,LogL,ConvError] = EM(tseries,settings,distvar)

% EM.m consists of the EM algorithm which is a tool for iterative maximum
% likelyhood estimation in order to estimate the variances that go into the
% matrices H (disturbance noise) and Q (process noise).
%
%
% HOW:     [var,LogL,ConvError] = EM(tseries,settings,var)
%
% Input:
%       tseries structure which should contain at least:
%               tseries.time [epochs x 1]           time vector
%               tseries.y [epochs x 1]              observation vector
%
%       settings, structure of which this function uses
%           settings.slope                  determines whether a slope will be estimated
%           settings.acc                    determines whether accelerations will be estimated
%           settings.cycle                  determines whether cycle components will be estimated
%           settings.numbercycles           number of cycle terms
%           settings.AR                     determines whether AR components will be estimated
%
%       distvar, structure with:
%           distvar.level                       variance of level term (usually zero to create smooth trend)
%           distvar.slope                       variance of slope term (zeta in Durbin&Koopman)
%           distvar.acc                         variance of acceleration
%           distvar.cycle [settings.numbercycles x 1] variance of cycle terms
%           distvar.irr                         variance of disturbance
%           distvar.AR                          variance of AR process
%
%       settings, structure of which this function uses
%           settings.maxiterEM              maximum number of iterations
%                                           (default = 1000)
%           settings.convEM                 convergence criterium, defined
%                                           as minimum difference between
%                                           subsequent log likelihood
%                                           (default = 1e-3)
%           settings.plotEM                 make plots
%
%       model, structure of which this function uses
%           model.W                         matrix containing the observed
%           dimensions per time step
%
%
% Output:
%       distvar, structure with:
%           distvar.level                       variance of level term (usually zero to create smooth trend)
%           distvar.slope                       variance of slope term (zeta in Durbin&Koopman)
%           distvar.acc                         variance of acceleration
%           distvar.cycle [settings.numbercycles x 1] variance of cycle terms
%           distvar.irr                         variance of disturbance
%           distvar.AR                          variance of A
%
%           LogL                            log likelihood
%
%           ConvError                       error messages
%
%
% note:
%
%
% Taco Broerse, Delft University of Technology, 2014
% d.b.t.broerse@tudelft.nl
%
%----------------------------------------------------------------------------
% Explanation of univariate structural time series model
%
% The theory behind the state space description, Kalman filtering
% and smoothing is largerly based on 'Time Series Analysis by State
% Space Methods' by J. Durbin and S.J. Koopman and
% to a minor extent on 'Forecasting, structural time series models and
% the Kalman filter' by A.C. Harvey
%
%----------------------------------------------------------------------------
% uses: StateSpace.m; Kalman.m; DisturbanceSmoother.m;
% InitialiseState.m; ExtractVarianceDisturbance.m;
%
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 June 2014 DBT Broerse
%
% Version 1.1 September 2014 DBT Broerse
% - added regression
% - added step intervention
%
% Version 1.2 October 2014 Thomas Frederikse
% - square root filter included
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
%
% Version 2.1 April 2015 DBT Broerse
% - no longer uses Q matrix but uses distvar.eta directly
%
% Version 2.2 November 2015 DBT Broerse / T Frederikse
% added AR process
%
% Version 2.3 June 2018 DBT Broerse
% added covariance for cycles
%
% Version 2.4 June 2018 DBT Broers
% extension to multivariate models
%
% Version 2.5 June 2019 DBT Broers
% extension to missing observations in multivariate models
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------
%
%%%% THEORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Formulation of the EM algorithm according to chapter 7.3.4 of Durbin and
% Koopman.
%
% The Expectation-Maximisation procedure consist of an iteration of
%
% VarEpsilon=VarEpsilon_tilde + 1/tseries.ntimes * VarEpsilon_tilde * sum_t=1,n
% (u_t^2 - D_t) * VarEpsilon_tilde
%
% where u_t and D_t are computed using the Kalman filter and Disturbance
% smoother
%
%
% VarEta = VarEta_tilde + 1/(tseries.ntimes-1) * VarEta_tilde *
% sum_t=1:n r_t-1*r_t-1' -N_t-1) * VarEta_tilde
%
% where r_t and N_t are computed using the Kalman filter and Disturbance
% smoother
%
% here Eta is the vector of process noise
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
narginchk(3,3);


%%

% defaults for maximum number of iterations
if ~any(strcmp('maxiterEM',fieldnames(settings)))
    settings.maxiterEM=1000;
elseif isempty(settings.maxiterEM)
    settings.maxiterEM=1000;
end
% defaults for convergence criterium (difference between subsequent
% loglikelihood values
if ~any(strcmp('convEM',fieldnames(settings)))
    settings.convEM=1e-3;
    disp('no settings.convEM')
elseif isempty(settings.convEM)
    settings.convEM=1e-3;
    disp('empty settings.convEM')
end
% defaults for plots
if ~any(strcmp('plotEM',fieldnames(settings)))
    settings.plotEM=true;
    disp('empty settings.plotEM')
end

LogL=[];
converged=false;
iteration=0;
ConvError.error = false;
%%
while (~converged)

    iteration=iteration+1;


    [model,distvar] = StateSpaceSetup(tseries,settings,distvar);

    % initialise
    if iteration==1
        % Setup matrices for state space formulation

        % set to full variance covariance matrix or irregular term
        VarEpsilon=distvar.varcovirr;

        % copy full variance covariance matrix
        VarEta=distvar.varcoveta;


    else
        % subsequent iterations

        % save structure
        [prevdistvar]=ExtractVarianceDisturbance(VarEta,VarEpsilon,model,tseries,settings);

        % set to full variance covariance matrix or irregular term
        VarEpsilon=distvar.varcovirr;

        % copy full variance covariance matrix, it removes covariances
        % between different components
        VarEta=distvar.varcoveta;

        % update disturbance variance matrix
        if settings.continuoustime
            for i=1:tseries.ntimes
                model.Q(:,:,i)=VarEta*tseries.dt(i);
            end
        else

            % integer time
            model.Q=VarEta;
        end

    end




    %% initialise state and variance

    [estfilter] = InitialiseStateIteration(tseries,settings,model,distvar);

    %% Run Kalman Filter (forward run, t=1:n)

    [estfilter] = Kalman(tseries,model,settings,estfilter);

    %% Disturbance smoothing (backward run, t=n:1)

    [estsmooth] = DisturbanceSmoother(tseries,model,settings,estfilter);


    %% Update variances of eta and epsilon

    % store previous values
    VarEpsilon_tilde = VarEpsilon;
    VarEta_tilde     = VarEta;

    % update of epsilon and eta, see Koopman - Disturbance smoother for state
    % space models (1993); equations 3.5, 3.6 and 3.7 (note the differerent
    % uses of G (our H) and H (our Q))

    if settings.multivariate
        % multivariate case
        SumTemp=zeros(tseries.ntseries);
        for i=1:tseries.ntimes
            SumTemp=SumTemp+estsmooth.u(:,i)*estsmooth.u(:,i)'-estsmooth.D(:,:,i);
        end
    else
        % univariate case
        SumTemp=zeros(tseries.ntseries);
        for i=1:tseries.ntimes
            SumTemp=SumTemp+estsmooth.u(i)*estsmooth.u(i)'-estsmooth.D(i);
        end
    end

    VarEpsilon = VarEpsilon_tilde + 1/tseries.ntimesobs * VarEpsilon_tilde * SumTemp * VarEpsilon_tilde;
    %%
    %VarEpsilon

    % multivariate case and univariate cases are equal
    SumTemp=zeros(model.dim);
    for i=2:tseries.ntimes
        SumTemp=SumTemp+estsmooth.r(:,i-1)*estsmooth.r(:,i-1)'-estsmooth.N(:,:,i-1);
    end


    VarEta = VarEta_tilde + 1/(tseries.ntimesobs-1) * VarEta_tilde * SumTemp * VarEta_tilde;

    %% compute log likelihood
    LogLPrev = LogL;
    [LogL]= LogLikelihood(tseries,estfilter,settings,model);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % store variances for plotting, in case of multivariate models, also
    % store covariances

    VarIrri{iteration}=distvar.irr;
    VarLeveli{iteration}=distvar.level;
    if settings.multivariate
        CovLeveli{iteration}=distvar.covlevel;
        CovIrri{iteration}=distvar.covirr;
    end
    LogLi(iteration)=LogL;
    if settings.slope
        VarSlopei{iteration}=distvar.slope;
        if settings.multivariate
            CovSlopei{iteration}=distvar.covslope;
        end
    end
    if settings.acc
        VarAcci{iteration}=distvar.acc;
        if settings.multivariate
            CovAcci{iteration}=distvar.covacc;
        end
    end
    if settings.cycle
        for i=1:settings.numbercycles
            VarCyclei{iteration,i}=distvar.cycle(:,i);
            VarCyclestari{iteration,i}=distvar.cyclestar(:,i);
            CovCyclestari{iteration,i}=distvar.covcyclestar(:,i);
        end
    end

    if settings.AR

        VarARi{iteration}=distvar.AR;

    end

    if settings.multivariate
        % also store covariances
        CovIrri{iteration}=distvar.covirr;
        CovLeveli{iteration}=distvar.covlevel;
        if settings.slope
            CovSlopei{iteration}=distvar.covslope;
        end
        if settings.acc
            CovAcci{iteration}=distvar.covacc;
        end
        % no covariances implemented for cycles yet

        if settings.AR
            CovARi{iteration}=distvar.covAR;
        end

    end


    %% check for convergence

    DiffLogL=(LogL-LogLPrev);

    
    
    if settings.MaxRelDeviationLogL~=0
        MaxDeviation=settings.MaxRelDeviationLogL*abs(LogL)*-1;
    else
        MaxDeviation=0;
    end


    if (DiffLogL<MaxDeviation) 
        warning('Likelihood decreases')
        fprintf('Using previous variance values  \n')
        fprintf(strcat('previous log likelihood: ',num2str(LogLPrev),'\n'))
        fprintf(strcat('current log likelihood:  ',num2str(LogL),'\n'))
        fprintf(strcat('iteration: ',num2str(iteration),'\n'))
        fprintf('If this happens in the beginning of the EM process, consider using different process noise variances \n')
        fprintf('Or reconsider a better model to describe the data \n')
        converged         = true;
        ConvError.error        = true;
        ConvError.type         = 'decreasing log likelihood';

        distvar=prevdistvar;
        ConvError.iter         = iteration;

    elseif (abs(DiffLogL)<settings.convEM)
        % convergence
        %  disp('solution converged')
        % extract individual variances from vectors, do not denormalise
        % state variable process noise

        [distvar]=ExtractVarianceDisturbance(VarEta,VarEpsilon,model,tseries,settings);



        ConvError.error        = false;
        ConvError.type         = 'log likelihood has converged';
        converged              = true;
        ConvError.iter         = iteration;
    else
        % no convergence yet, proceed

        % put updated VarEta & VarEpsilon back in distvar
        [distvar]=ExtractVarianceDisturbance(VarEta,VarEpsilon,model,tseries,settings);


    end



    % check whether maximum nr of iterations has been reached
    if iteration==settings.maxiterEM
        disp(strcat('reached ',num2str(settings.maxiterEM),' iterations, stopping now'))
        converged=true;
        ConvError.error='true';
        ConvError.type='no convergence';
        ConvError.iter         = iteration;

    end

end



%% plot convergence of estimated parameters
if settings.plotEM
    MarkerSize=10;
    LineWidth=1.5;
    NumSubPanels=3;
    if settings.slope
        NumSubPanels=NumSubPanels+1;
    end
    if settings.acc
        NumSubPanels=NumSubPanels+1;
    end
    if settings.cycle
        for i=1:settings.numbercycles
            NumSubPanels=NumSubPanels+1;
        end
    end
    if settings.AR
        NumSubPanels=NumSubPanels+1;
    end


    % determine number of sub panels
    nCols=ceil(sqrt(NumSubPanels));
    nRows=ceil(NumSubPanels/nCols);

    fig=figure;hold on
    fig.Position=[100 100 800 800];
    % get colors
    ColorOrder=get(gca,'ColorOrder');
    nColors=length(ColorOrder);

    for p=1:tseries.ntseries*(model.dim+1)
        Colors(p,:)=ColorOrder(mod(p-1,nColors)+1,:);
    end

    %% log likelihood
    iPlot=1;
    subplot(nRows,nCols,iPlot)
    semilogy([1:iteration],LogLi,'-','LineWidth',LineWidth);hold on
    xlim([0 iteration])
    ylim([prctile(LogLi,1) max(LogLi)])
    xlabel('iteration')
    title('log likelihood')


    iPlot=iPlot+1;
    %% level
    subplot(nRows,nCols,iPlot)


    for p=1:tseries.ntseries
        levelvec=zeros(iteration,1);
        for i=1:iteration
            levelvec(i)=VarLeveli{i}(p);
        end
        if VarLeveli{i}(end) ~=0
            semilogy([1:iteration],levelvec,'-','MarkerSize',MarkerSize,'Color',Colors(p,:),'LineWidth',LineWidth);hold on
        else
            plot([1:iteration],levelvec,'-','Color',Colors(p,:),'LineWidth',LineWidth); hold on
        end

    end
    xlim([0 iteration])

    title('level variance')
    %% irregular
    iPlot=iPlot+1;
    subplot(nRows,nCols,iPlot)


    ilegend=0;
    for p=1:tseries.ntseries
        irrvec=zeros(iteration,1);
        for i=1:iteration
            irrvec(i)=VarIrri{i}(p);
        end

        ilegend=ilegend+1;
        h(ilegend)=semilogy([1:iteration],irrvec,'-','MarkerSize',MarkerSize,'Color',Colors(p,:),'LineWidth',LineWidth); hold on
        legendstr{ilegend}=strcat('var serie:',num2str(p));


        if settings.multivariate
            for pp=p+1:tseries.ntseries
                corrirrvec=zeros(iteration,1);
                for i=1:iteration
                    corrirrvec(i)=CovIrri{i}(p,pp)/sqrt(VarIrri{i}(p)*VarIrri{i}(pp));
                end
                ilegend=ilegend+1;
                yyaxis right
                h(ilegend)=plot([1:iteration],corrirrvec,'--','Color',Colors(tseries.ntseries+pp,:),'LineWidth',LineWidth);
                axesright=gca;
                axesright.YColor=[0 0 0];
                ylim([-1 1])
                ylabel('correlation')
                yyaxis left
                legendstr{ilegend}=strcat('corr serie:',num2str(p),':',num2str(pp));
            end
        end

    end
    xlim([0 iteration])
    title('irregular variance')
    legend(h,legendstr,'Location','Best')
    legend('boxoff')



    if settings.slope
        iPlot=iPlot+1;
        subplot(nRows,nCols,iPlot)


        for p=1:tseries.ntseries
            slopevec=zeros(iteration,1);
            for i=1:iteration
                slopevec(i)=VarSlopei{i}(p);
            end
            semilogy([1:iteration],slopevec,'-','MarkerSize',MarkerSize,'Color',Colors(p,:),'LineWidth',LineWidth)
            hold on
            if settings.multivariate
                for pp=p+1:tseries.ntseries
                    corrslopevec=zeros(iteration,1);
                    for i=1:iteration
                        corrslopevec(i)=CovSlopei{i}(p,pp)/sqrt(VarSlopei{i}(p)*VarSlopei{i}(pp));
                    end

                    yyaxis right
                    plot([1:iteration],corrslopevec,'--','Color',Colors(tseries.ntseries+pp,:),'LineWidth',LineWidth);
                    ylim([-1.05 1.05])
                    ylabel('correlation')
                    axesright=gca;
                    axesright.YColor=[0 0 0];
                    yyaxis left

                end

            end
        end
        xlim([0 iteration])
        title('slope variance')
    end

    %% acceleration

    if settings.acc
        iPlot=iPlot+1;
        subplot(nRows,nCols,iPlot)


        for p=1:tseries.ntseries
            accvec=zeros(iteration,1);
            for i=1:iteration
                accvec(i)=VarAcci{i}(p);
            end
            semilogy([1:iteration],accvec,'-','MarkerSize',MarkerSize,'Color',Colors(p,:),'LineWidth',LineWidth); hold on

        end
        xlim([0 iteration])
        title('acceleration variance')

    end

    %% cycles
    if settings.cycle

        for i=1:settings.numbercycles
            iPlot=iPlot+1;
            subplot(nRows,nCols,iPlot)
            ilegend=0;
            legendstr=[];
            clear h
            for p=1:tseries.ntseries
                if VarCyclei{end}(p)~=0
                    cyclevec=zeros(iteration,1);
                    cyclestarvec=zeros(iteration,1);
                    cyclecorvec=zeros(iteration,1);
                    for ii=1:iteration
                        cyclevec(ii)=VarCyclei{ii,i}(p);
                        cyclestarvec(ii)=VarCyclestari{ii,i}(p);
                        cyclecorvec(ii)=CovCyclestari{ii,i}(p)/sqrt(cyclevec(ii)*cyclestarvec(ii));
                    end
                    ilegend=ilegend+1;
                    h(ilegend)=semilogy([1:iteration],cyclevec,'-','MarkerSize',MarkerSize,'Color',Colors(p,:),'LineWidth',LineWidth); hold on
                    legendstr{ilegend}=strcat('var serie:',num2str(p));
                    semilogy([1:iteration],cyclestarvec,'-','Color',Colors(p,:),'LineWidth',LineWidth)
                    yyaxis right

                    ilegend=ilegend+1;
                    h(ilegend)=semilogy([1:iteration],cyclecorvec,'--','Color',Colors(tseries.ntseries+p,:),'LineWidth',LineWidth);
                    legendstr{ilegend}=strcat('corr c c^* serie:',num2str(p));
                    ylim([-1.05 1.05])
                    ylabel('correlation')
                    axesright=gca;
                    axesright.YColor=[0 0 0];
                    yyaxis left
                else
                    cyclevec=zeros(iteration,1);
                    cyclestarvec=zeros(iteration,1);
                    cyclecorvec=zeros(iteration,1);
                    for ii=1:iteration
                        cyclevec(ii)=VarCyclei{ii,i}(p);
                        cyclestarvec(ii)=VarCyclestari{ii,i}(p);
                        cyclecorvec(ii)=CovCyclestari{ii,i}(p)/sqrt(cyclevec(ii)*cyclestarvec(ii));
                    end
                    ilegend=ilegend+1;
                    h(ilegend)=plot([1:iteration],cyclevec,'-','MarkerSize',MarkerSize,'Color',Colors(p,:),'LineWidth',LineWidth); hold on
                    legendstr{ilegend}=strcat('var serie:',num2str(p));
                    plot([1:iteration],cyclestarvec,'.','Color',Colors(p,:),'LineWidth',LineWidth)
                    yyaxis right
                    ilegend=ilegend+1;
                    h(ilegend)=plot([1:iteration],cyclecorvec,':','Color',Colors(p,:),'LineWidth',LineWidth);
                    legendstr{ilegend}=strcat('corr c c^* serie:',num2str(p));
                    ylim([-1.05 1.05])
                    ylabel('correlation')
                    axesright=gca;
                    axesright.YColor=[0 0 0];
                    yyaxis left
                end
            end
            legend(h,legendstr,'Location','Best')
            legend('boxoff')
            xlim([0 iteration])
            title(strcat('(co)var cycle, period:',num2str(settings.periodscycle(i))))
            %         end
        end
    end

    %% AR
    if settings.AR
        iPlot=iPlot+1;
        subplot(nRows,nCols,iPlot)
        for p=1:tseries.ntseries
            ARvec=zeros(iteration,1);
            for i=1:iteration
                ARvec(i)=VarARi{i}(p);

            end
            semilogy([1:iteration],ARvec,'-','MarkerSize',MarkerSize,'Color',Colors(p,:),'LineWidth',LineWidth); hold on
            if settings.multivariate
                for pp=p+1:tseries.ntseries
                    corrARvec=zeros(iteration,1);
                    for i=1:iteration
                        corrARvec(i)=CovARi{i}(p,pp)/sqrt(VarARi{i}(p)*VarARi{i}(pp));
                    end

                    yyaxis right
                    plot([1:iteration],corrARvec,'--','Color',Colors(tseries.ntseries+pp,:));
                    ylim([-1.05 1.05])
                    ylabel('correlation')
                    axesright=gca;
                    axesright.YColor=[0 0 0];
                    yyaxis left

                end

            end


        end
        xlim([0 iteration])
        title(strcat('AR(1) variance for phi:',num2str(settings.ARphi)))
    end

    % add title
    if isfield(tseries,'name')
        if ConvError.error
            convstr=ConvError.type;
        else
            convstr='converged';
        end
        titlestr=strcat({'EM variance optimization for '},{tseries.generalname},{' '},{strcat('result:',convstr)});
        annotation('textbox', [0 0.9 1 0.1], ...
            'String', titlestr, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center')

    end

    % save figure
    if settings.saveplots
        if isfield(tseries,'name')
            if isfield(settings,'maindir')
            filename=strcat(strcat(settings.maindir,'/',settings.savedir,'/','EM_optimization_',tseries.generalname));
            else
            filename=strcat(strcat(settings.savedir,'/','EM_optimization_',tseries.generalname));

            end
        else
            filename=strcat(strcat(settings.maindir,'/',settings.savedir,'/','EM_optimization'));
        end
        savefig(filename);

        if settings.savepng

            print(fig,strcat(filename,'.png'),'-dpng')
        end
    end



end



