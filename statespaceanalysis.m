function[tseries,model,distvar,settings,estfilter,fitfilter,estsmooth,fitsmooth,varsmooth,ConvError,covsmooth,meanslope]=statespaceanalysis(tseries,settings,distvar);

% statespaceanalysis runs complete state space analysis for time series
% tseries
%
% HOW:     [tseries,model,distvar,settings,estfilter,fitfilter,estsmooth,fitsmooth,varsmooth,ConvError]=statespaceanalysis(tseries,settings,distvar);


% Input:
%           tseries  times series structure, which should contain:
%               tseries.time [epochs x 1]           time vector
%               tseries.y [epochs x 1]              observation vector
%
%           settings                    program and model settings structure
%           distvar                     disturbance variances of process and irregular component
%
% Output:
%
%           tseries                     time series structure
%           model                       state space model
%           distvar                     disturbance variances (updated)
%           settings                    program and model settings structure
%           estfilter                   filter output (including state
%                                       vector and error variance)
%           fitfilter                   reconstructed fits from filter
%           estsmooth                   smoother output (including state
%                                       vector and error variance)
%           fitsmooth                   reconstructed fits from smoother
%           varsmooth                   error variances state vector of smoother
%           ConvError                   error information from parameter estimation
%           covsmooth                   smoothed state covariances
%           meanslope                   mean slope
%
%
%
%%%%% Function to run time series analysis based on state space
%%%%% formulation and makes use of Kalman filtering.
%%%%% Purpose is to estimate time variable seasonal (cyclic in state space
%%%%% nomenclature) signals and time variable trends
%%%%%
%%%%% The program consists of a few essential discrete steps:
%%%%% 0. User selection of the state space discription of the model that has
%%%%%    to be fit to the observations and charachteristics of disturbance
%%%%%    and process variances.
%%%%% 1. Estimation of irregular (unmodeled signal) and process (model)
%%%%%    variances using maximalisation of the log likelihood (this
%%%%%    includes an iteration of steps 1,2,3,5.
%%%%% 2. Creatinging matrices for later use in Kalman filter: design matrix
%%%%%    Z, transition matrix T, process noise Q and observational noise H
%%%%% 3. Running Kalman filter, which carries out a recursive update of the
%%%%%    state (forward loop)
%%%%% 4. State smoothing (backward loop), which updates the state using all
%%%%%    available observations
%%%%% 5. Disturbance smoothing (backward loop), which computes smoothed
%%%%%    estimates of the disturbance and process noise using all available
%%%%%    observations.
%%%%%
%%%%% Furthermore, this script makes some plots of results and inputs
%%%%%
%%%%% The theory behind the state space description and Kalman filtering
%%%%% and smoothing is largerly based on 'Time Series Analysis by State
%%%%% Space Methods' by J. Durbin and S.J. Koopman and
%%%%% to a minor extent on 'Forecasting, structural time series models and
%%%%% the Kalman filter' by A.C. Harvey in case lag operators are chosen
%%%%% for the state space description.
%%%%%
%%%%%
%%%%%
%
% Taco Broerse, Delft University of Technology, 2014
% d.b.t.broerse@tudelft.nl
%
%
%--------------------------------------------------------------------------
% Requires the following additional files:
%
% StateSpacePreProcess.m    preprocessing (normalising) of time vector
% EM.m                      Estimation/Maximisation algorithm
% StateSpaceSetup.m              defines the time series model in the state space frame
% Kalman.m                  runs kalman filter (forward pass)
% Smoother.m                smoothes time series (backward pass)
% DisturbanceSmoother.m     smoothes disturbance and process noise (backward pass)
% PostProcessState.m        post-processing of state
% ExtractVariance.m         subroutine to extract process variances from vector
% LogLikeliHood.m           computes log likelihood
% InitialiseState.m         initialises the state for the first epoch (no
%                           observations included)
% FindMissingEpochs.m       find missing epochs in homogenously spaced data
%                           and insert these epochs with NaN in the
%                           observation vector
%
%
%--------------------------------------------------------------------------
% revision history
%
% Version 1.0 October 2014 DBT Broerse
%
% Version 1.1 June 2018 DBT Broerse, now includes AR process as implemented
% by Thomas Frederikse
%
% Version 1.2 June 2018 DBT Broerse, includes check for disturbance
% variances
%
% Version 1.3 May 2020 DBT Broerse, includes optimization of different AR phi
% parameters for each time series
%
% Version 1.4 bug in missing epoch detection corrected
%
%--------------------------------------------------------------------------
%%%%% Still to be implemented:
%%%%%
%%%%% - diffuse initialisation
%--------------------------------------------------------------------------

%defaults
ConvError=[];

% make directory for saving
if isfield(settings,'maindir')
    savedir=(strcat(settings.maindir,'/',settings.savedir));
else
    savedir=(settings.savedir);
end
if ~isdir(savedir)
    mkdir(savedir)
    
end


%% find missing epochs
% checkmissingepoch=1;
if settings.continuoustime == false && isempty(tseries.missingepochs)
    %     checkmissingepoch=1;
    % end
    % if settings.multivariate
    %     % also check missing epochs for multivariate series
    %     checkmissingepoch=1;
    % end
    % if checkmissingepoch
    tseries.oldtime=tseries.time;
    tseries.oldy=tseries.y;
    [tseries]=FindMissingEpochs(tseries);
    if ~ isempty(tseries.missingepochs)
        if settings.checkplots
            figure;hold on
            plot(tseries.oldtime,tseries.oldy)
            plot(tseries.time,tseries.y,'o')
            title('check for inserting missing observations with NaN')
            
            vline(timecomplete(tseries.missingepochs))
        end
    else
        clear tseries.oldtime tseries.oldy
    end
else
    if settings.multivariate
        ntseries=numel(tseries.y)/numel(tseries.time);
        for p=1:ntseries
            tseries.missingepochs{p}=find(isnan(tseries.y(:,p)));
        end
    else
        if ~isempty(find(isnan(tseries.y)))
            disp('note that there are missing epochs in a univariate series')
            tseries.missingepochs=find(isnan(tseries.y));
            
        end
    end
end



%% Normalize time

[tseries,settings]=StateSpacePreProcess(tseries,settings);

%% check
if ~isfield(tseries,'generalname')
    if tseries.ntseries==1
        tseries.generalname=tseries.name;
    else
        tseries.generalname=cell2mat(tseries.name);
    end
end


%% remove epochs that are simultaneous with steps (interventions)

if settings.intervention
    
    if settings.removeinterventionepochs
        tseries.yremoved=[];tseries.timeremoved=[];
        irem=0;
        % loop on steps (interventions)
        for i=1:length(settings.interventiontime)
            
            % check the time difference with the intervention time
            difftime=abs(tseries.time-settings.interventiontime(i));
            settings.removeperiod=0.9/365; % ie remove everything within that day
            
            % find indices of epochs that should be removed
            indexremove=find(difftime < settings.removeperiod);
            disp(['remove indices:' indexremove])
            
            if ~isempty(indexremove)
                for ii=length(indexremove):-1:1
                    irem=irem+1;
                    
                    %       disp(strcat('removing time step:',num2str(tseries.time(indexremove(ii)))))
                    tseries.yremoved(:,irem)=tseries.y(:,indexremove(ii));
                    tseries.timeremoved(irem)=tseries.time(indexremove(ii));
                    % remove epochs
                    tseries.time(indexremove(ii))=[];
                    tseries.y(:,indexremove(ii))=[];
                end
                
            end
            
        end
        figure; hold on
        title('removed epochs (in red)')
        npanels=tseries.ntseries;
        nrows=ceil(sqrt(npanels));
        ncols=ceil(npanels/nrows);
        MarkerSize=10;
        for i=1:tseries.ntseries
            subplot(nrows,ncols,i); hold on
            plot(tseries.time,tseries.y(i,:),'.k','MarkerSize',MarkerSize)
            plot(tseries.timeremoved,tseries.yremoved(i,:),'r.','MarkerSize',MarkerSize*1.5)
        end
        
        
        [tseries,settings]=StateSpacePreProcess(tseries,settings);
    end
end


%% Do checks on disturbance variances and if AR is present, length of AR phi coefficient vector

[distvar,settings]=ProcessNoiseCheck(tseries,settings,distvar);

%% START ACTUAL STATE SPACE ANALYSIS
% Estimate process variance and irregular component variance

if (~settings.fixprocvar)
    if strcmp(settings.optvar,'EM')
        if settings.AR&&settings.optAR
            if settings.ARorder>1 %% DOES NOT WORK FOR AR(i>1)!!!
                error('optimization of phi not possible for AR(i>1)')
            end
            
            % in order to prevent a high number of plots
            originalplotsetting=settings.plotEM;
            settings.plotEM=0;
            
            
            % Do global optimization of the loglikelihood for phi [0 1] with EM
            % algorithm to determine disturbance variances.
            % loop on a range of phi values, and take the one that gives the
            % best log likelihood
            if settings.ARphiMulti
                if tseries.ntseries > 2
                    error('optimization of phi separately, not implemented for more than two time series')
                end
                % currently only works for two time series
                lengthARphiRange=length(settings.ARphiRange);
                logl=zeros(lengthARphiRange,lengthARphiRange);
            else
                lengthARphiRange=length(settings.ARphiRange);
                logl  = zeros(lengthARphiRange,1);
            end
            
            disp('maximizing log likelihood for various values of phi (AR(1))')
            
            if settings.ARphiMulti && settings.multivariate
                % case: different phi for all time series
                for i=1:numel(settings.ARphiRange)
                    for j=1:numel(settings.ARphiRange)
                        disp(strcat('phi 1:',num2str(settings.ARphiRange(i)),' phi 2:',num2str(settings.ARphiRange(j))))
                        % reset to original initial values
                        
                        settings.ARphi = [settings.ARphiRange(i) ; settings.ARphiRange(j)];
                        [distvartemp,LogL,~] = EM(tseries,settings,distvar);
%                         if abs(mean(distvartemp.corrcyclestar)-1)> eps
%                             % do this check to see if the cycle terms stay
%                             % correlated
%                         distvartemp.corrcyclestar
%                         end
                        logl(i,j)=LogL;
                        disp(strcat('var slope:',num2str(distvartemp.slope),' var AR:',num2str(distvartemp.AR)))
                        disp(strcat('log likelihood:',num2str(LogL)))
                        disp(' ')
                    end
                end
                if tseries.ntseries > 2
                    error('not implemented yet for more than 2 time series')
                end
                % take the value pertaining to the maximum log likelihood
                [a,b]=find(logl==max(logl(:)));
                settings.ARphi=[];
                settings.ARphi(1,1)=settings.ARphiRange(a);
                settings.ARphi(2,1)=settings.ARphiRange(b);
                
            else
                % case: same phi for all time series
                for i=1:numel(settings.ARphiRange)
                    disp(strcat('phi:',num2str(settings.ARphiRange(i))))
                    % reset to original initial values
                    
                    settings.ARphi = settings.ARphiRange(i);
                    [distvartemp,LogL,~] = EM(tseries,settings,distvar);
                    
                    logl(i)=LogL;
                    disp(strcat('var slope:',num2str(distvartemp.slope),' var AR:',num2str(distvartemp.AR)))
                end
                
                % take the value pertaining to the maximum log likelihood
                settings.ARphi = settings.ARphiRange(logl==max(logl(:)));
            end
            
            
            % reset to original setting
            settings.plotEM=originalplotsetting;
            
            if settings.plotEM
                
                
                fig=figure; hold on
                if settings.ARphiMulti && settings.multivariate
                    % don't have a method yet to show more than 2
                    % dimensions
                    colormaps=load('lajolla.mat');
                    cmap=flipud(colormaps.lajolla);
                    
                    imagesc(settings.ARphiRange,settings.ARphiRange,(logl))
                    plot(settings.ARphi(2),settings.ARphi(1),'ko')
                    xlabel('\phi dimension 2')
                    ylabel('\phi dimension 1')
                    c=colorbar;
                    c.Label.String='log likelihood';
                    axis tight
                    axis equal
                    colormap(cmap)
                else
                    plot(settings.ARphiRange,logl,'.')
                    plot(settings.ARphi,max(logl),'o')
                    xlabel('\phi AR')
                    ylabel('log likelihood')
                end
                
                title('optimization of phi')
                
                % save figure
                if settings.saveplots
                    if isfield(tseries,'name')
                        filename=strcat(settings.maindir,'/',settings.savedir,'/','AR_phi_optimization_',tseries.generalname);
                    else
                        filename=strcat(settings.maindir,'/',settings.savedir,'/','AR_phi_optimization');
                    end
                    savefig(filename);
                    if settings.savepng
                        print(fig,strcat(filename,'.png'),'-dpng')
                    end
                end
            end
            
        end
        
        % run EM algorithm for optimizing process noise and irregular component
        % variance
        [distvar,logl,ConvError] = EM(tseries,settings,distvar);
        
        if settings.optAR
            % save output of testing AR phi
            ConvError.LogLfuncARphi=logl;
        end
    else
        settings.optvar
        error('only optimization algorithm implemented is EM')
    end
end


%% Setup matrices for state space formulation

[model,distvar] = StateSpaceSetup(tseries,settings,distvar);


%% Initialise state vector and variance

[estfilter] = InitialiseStateIteration(tseries,settings,model,distvar);

%% Run Kalman Filter (forward run, t=1:n)

[estfilter] = Kalman(tseries,model,settings,estfilter);
%% extract fit
[fitfilter] = PostProcessState(estfilter,settings,model,tseries);

%% State Smoothing (backward run, t=n:1)

[estsmooth] = Smoother(tseries,settings,model,estfilter);

%% Disturbance smoothing (backward run, t=n:1)

[estsmooth] = DisturbanceSmoother(tseries,model,settings,estfilter,estsmooth);

%% log likelihood
[LogL]= LogLikelihood(tseries,estfilter,settings,model);
estfilter.LogL=LogL;

%% Post-process
% get time series of components and cycle information (phase and amplitude)
[fitsmooth] = PostProcessState(estsmooth,settings,model,tseries);% get time series of state error variances

%% Post-process variances
% always denormalise
denormalise=true;
[varsmooth]=ExtractVarianceState(estsmooth.V,tseries,settings,model,denormalise);

%% Covariances

if settings.estimatemeanslope
    [covsmooth] = Covariances(tseries,estfilter,estsmooth,model);
else
    covsmooth = [];
end
%% Mean trend
if settings.estimatemeanslope
    [meanslope] = PostProcessTrend(settings,tseries,fitsmooth,varsmooth,covsmooth,model);
else
    meanslope = [];
end

end
