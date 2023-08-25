%%%%% Example script to run Statespace. This program performs a time series
%%%%% analysis based on state space formulation and makes use of Kalman filtering.
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
%%%%% Still to be implemented:
%%%%%
%%%%% - diffuse initialisation
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
% LeastSquaresTrend         deterministic trend analysis using least squares
% SpectralAnalysis          frequency analysis of observation and estimates
% EM.m                      Estimation/Maximisation algorithm
% StateSpace.m              defines the time series model in the state space frame
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
%                           observation vector.
%
% for time series with varying time step a few different files are used instead:
%
% lomb.m                    Lomb-Scargle periodigram
%
%
%--------------------------------------------------------------------------
% revision history
%
% Version 1.0 June 2014 DBT Broerse
%
% Version 2.0 October 2014 DBT Broerse
% - improved I/O
% - state space analysis now in subroutine statespace
%
% Version 3.0 June 2018 DBT Broerse
% include multivariate models
%
% Version 3.1 March 2022 DBT Broerse
% inclusion multiple data test files and extended statistical testing
%--------------------------------------------------------------------------

close all
clear all;   clc


%%
restoredefaultpath
% directory for running state space, change to your own directory
maindir = '/Users/tacobroerse/MATLAB/work/state_space/working_folder_multivariate/test/';
cd(maindir)
% directory for state space code, change to your own directory
codedir='/Users/tacobroerse/MATLAB/work/state_space/working_folder_multivariate/';

% update path
addtopath={codedir,strcat(codedir,'/PlotFunctions/'),strcat(codedir,'/colormaps/')};



for i=1:length(addtopath)
    path(path,addtopath{i})
end


%% set case (to test different time series)

%DataType = 'GRACE';
DataType = 'GRACE-w-accelerations';  % same, but with stochastic slope replaced by stochastic acceleration

%DataType = 'Nile';

% DataType = 'TideGauge';
%
%DataType = 'synthetic';
%
% DataType = 'InSAR'

% DataType = 'InSARHarlingen1'
%
% DataType = 'InSARHarlingen2'

%DataType = 'KSI_UK';

%DataType = 'Paper_Willen_etal_2021_basin_21';

%% settings

% structure settings contains options and settings that do not change
% defaults
[settings,distvar]=SetDefaults;

settings.multivariate = 0;% set 0 for univariate time series, 1 for multivariate time series
settings.savedir=strcat('results/',DataType,'/');% folder for saving results
settings.saveplots=1;
settings.savepng=1;
settings.batchname=DataType;
settings.maindir=maindir;
settings.saveresults=1;

settings.estimatemeanslope=1;% compute mean slope
settings.type = 'Normal_Slope_and_Acc'; % type of state space formulation ('Normal_Slope_and_Acc' is default, you may use StateSpaceType=[])
settings.filter = 'DurbinKoopman'; % kalman filter and smoother type (default is 'DurbinKoopman', other option SquareRoot)
settings.checktime = false; % check whether time is integer
settings.fixprocvar = false; % if true no variances for the irregular component and process noise are estimated
settings.checkplots = false; % produce extra plots for checks
settings.lsq = 1; % do a least squares analysis for comparison
settings.maxiterEM = 1000; % Iteration maximum of EM algorithm
settings.convEM = 1e-3; % Convergence criterium EM algorithm (min difference in log likelihood between two consequtive iterations)
settings.convEM = 1e-4;
settings.optvar = 'EM'; % variance optimizer
settings.normalise= true; % do or do not normalise time (false only possible with continuous time)
settings.largeP = 1e4;% large initial variance
settings.distvarcovariancecycles=0;

%% Data should be provided in columns, first column time, second column data
% Also, set which components should be included in the state space model (and thus estimated)

if strcmp(DataType,'GRACE') || strcmp(DataType,'GRACE-w-accelerations')
    DataFile = 'data/china_trend.txt';
    settings.slope = true; %slope, for a local linear trend model a.k.a. integrated random walk if stochastic (process variance>0), linear drift if deterministic (process variance=0), keep at true for smooth trends
    settings.acc = false; %acceleration, however seems not so useful when trend is stochastic
    if strcmp(DataType,'GRACE-w-accelerations')
        settings.acc = true;
    end
    settings.acclsq = false; % acceleration for least squares
    settings.cycle = true; %cycle terms (harmonic series)
    settings.regression = false; % regression (fitting functions)
    settings.intervention = false; % step function at time tau
    settings.continuoustime = true; % if false time steps are equal
    
    % Autoregressive term
    settings.AR      = false; % Use autoregressive component
    settings.ARorder = 1;    % Choose order of AR-component
    settings.ARphi   = [0.1772]; % Set autoregressive coefficients
    settings.optAR   = true; % Optimize choice of autoregressive component. Use only for order-1 AR processes
    settings.ARphiRange = [0:0.1:1];
elseif strcmp(DataType,'TideGauge')
    DataFile = 'data/thomas_tidegauge2.mat';
    RegressorFile = 'data/regressors_thomas.mat';
    settings.slope = true; %slope, for a local linear trend model a.k.a. integrated random walk if stochastic (process variance>0), linear drift if deterministic (process variance=0), keep at true for smooth trends
    settings.acc = false; %acceleration, however seems not so useful when trend is stochastic
    settings.acclsq = false; % acceleration for least squares
    settings.cycle = true; %cycle terms (harmonic series)
    settings.regression = true; % regression (fitting functions)
    settings.intervention = false; % step function at time tau
    settings.continuoustime = true; % if false time steps are equal
    
    % Autoregressive term
    settings.AR      = true; % Use autoregressive component
    settings.ARorder = 1;    % Choose order of AR-component
    settings.ARphi   = [0.1772]; % Set autoregressive coefficients
    settings.optAR   = true; % Optimize choice of autoregressive component. Use only for order-1 AR processes
    settings.ARphiRange = [0:0.1:1];
elseif strcmp(DataType,'Nile')
    DataFile = 'data/nile.txt';
    settings.slope = false; %slope, for a local linear trend model a.k.a.
    %                       integrated random walk if stochastic (process variance>0),
    %                       linear drift if deterministic (process variance=0),
    %                       keep at true for smooth trends; false turns off
    %                       the slope
    settings.acc = false; %acceleration, however seems not so useful when trend is stochastic
    settings.acclsq = false; % acceleration for least squares
    settings.cycle = false; %cycle terms (harmonic series)
    settings.regression = false; % regression (fitting functions)
    settings.intervention = false; % step function at time tau
    settings.continuoustime = false; % if false time steps are equal
elseif strcmp(DataType,'InSAR')
    DataFile = 'data/InSAR_max_cycle.mat';
    settings.slope = true; %slope, for a local linear trend model a.k.a. integrated random walk if stochastic (process variance>0), linear drift if deterministic (process variance=0), keep at true for smooth trends
    settings.acc = false; %acceleration, however seems not so useful when trend is stochastic
    settings.acclsq = false; % acceleration for least squares
    settings.cycle = true; %cycle terms (harmonic series)
    settings.regression = false; % regression (fitting functions)
    settings.intervention = false; % step function at time tau
    settings.continuoustime = true; % if false time steps are equal
elseif strcmp(DataType,'InSARHarlingen1')
    DataFile = 'data/Harlingendecrease.mat';
    settings.slope = true; %slope, for a local linear trend model a.k.a. integrated random walk if stochastic (process variance>0), linear drift if deterministic (process variance=0), keep at true for smooth trends
    settings.acc = false; %acceleration, however seems not so useful when trend is stochastic
    settings.acclsq = true; % acceleration for least squares
    settings.cycle = false; %cycle terms (harmonic series)
    settings.regression = false; % regression (fitting functions)
    settings.intervention = false; % step function at time tau
    settings.continuoustime = true; % if false time steps are equal
elseif strcmp(DataType,'InSARHarlingen2')
    DataFile = 'data/Harlingenincrease.mat';
    settings.slope = true; %slope, for a local linear trend model a.k.a. integrated random walk if stochastic (process variance>0), linear drift if deterministic (process variance=0), keep at true for smooth trends
    settings.acc = false; %acceleration, however seems not so useful when trend is stochastic
    settings.acclsq = true; % acceleration for least squares
    settings.cycle = false; %cycle terms (harmonic series)
    settings.regression = false; % regression (fitting functions)
    settings.intervention = false; % step function at time tau
    settings.continuoustime = true; % if false time steps are equal
elseif strcmp(DataType,'synthetic')
    settings.slope = true; %slope, for a local linear trend model a.k.a. integrated random walk if stochastic (process variance>0), linear drift if deterministic (process variance=0), keep at true for smooth trends
    settings.acc = false; %acceleration, however seems not so useful when trend is stochastic
    settings.acclsq = false; % acceleration for least squares
    settings.cycle = true; %cycle terms (harmonic series)
    settings.regression = false; % regression (fitting functions)
    settings.intervention = false; % step function at time tau
    settings.continuoustime = false; % if false time steps are equal
elseif strcmp(DataType,'KSI_UK')
    DataFile = 'data/KSI_UK.mat';
    settings.slope = false;
    settings.acc = false;
    settings.acclsq = false;
    settings.cycle = true;
    settings.regression = false;
    settings.intervention = false;
    settings.continuoustime = true;
elseif strcmp(DataType,'Paper_Willen_etal_2021_basin_21')
    DataFile = 'data/Willen_etal_basin_21.txt';
    settings.lsq=0;
    settings.slope = true;
    settings.acc = false;
    settings.acclsq = false;
    settings.cycle = true;
    settings.regression = false;
    settings.intervention = false;
    settings.continuoustime = true;
    settings.multivariate = true;
    settings.startARat0=true;% start AR process at 0
    settings.distvarcovariancecycles=0; % no covariance for the cycles
    settings.ARphiMulti=1;% separate phi value for each time series
    settings.multivariatedistcovarianceslope=0;% use covariances in the slope disturbance variance
     % Autoregressive term
    settings.AR      = true; % Use autoregressive component
    settings.ARorder = 1;    % Choose order of AR-component
    settings.ARphi   = [0.8]; % Set autoregressive coefficients
    settings.optAR   = true; % Optimize choice of autoregressive component. Use only for order-1 AR processes
    %settings.optAR=false;
    settings.ARphiRange = [0.5:0.1:1];
end


% periods of cycles, i.e. annual: 1, semi-annual: 0.5, etc. in case time is
% in years
if strcmp(DataType,'GRACE') || strcmp(DataType,'GRACE-w-accelerations')
    settings.periodscycle = [1 0.5];
elseif strcmp(DataType,'TideGauge')
    settings.periodscycle = [1 0.5];
elseif strcmp(DataType,'KSI_UK')
    settings.periodscycle = 1;
elseif strcmp(DataType,'Paper_Willen_etal_2021_basin_21')
    settings.periodscycle = [1 0.5];
else
    settings.periodscycle = [1];
end



if settings.intervention
    if strcmp(DataType,'synthetic')
        settings.interventiontime=10;
    else
        settings.interventiontime=[2007.6961 ];
        % settings.interventiontimes=[2.003208333333333e+03+1];
    end
else
    settings.interventiontime=[];
end

% set regressors to empty if no regression is needed
if ~settings.regression
    tseries.regressors=[];
end

tseries.missingepochs= []; % supply indexes of time vector with missing observations (optional and only for equal time steps)

%% Set (initial) process variances for each component
% Process variances determine the amount of change per unit time step in
% each component and can be adjusted or estimated.
%
% Setting all process variances to zero turns estimation into a deterministic one.
% Also a mixed approach with deterministic and stochastic components is
% possible.
%
% Value of variance should be adjusted such that ideally each component (cyclic,
% trend) absorbs the right frequencies of the original signal and the
% residuals only contain high frequency noise. This can be chosen manually
% or estimated by maximizing the log likelihood. When estimated the values
% here will serve as initial values in an iterative approach. Variances set
% to zero will stay zero when variances are estimated.

% structure var contains variance values that may change
clear distvar
if strcmp(DataType,'GRACE') || strcmp(DataType,'GRACE-w-accelerations')
    if ~settings.normalise
        scale=1/(0.0833);
    else
        scale=1;
    end
    %deviation=0.95;
    distvar.irr=4*scale; % variance for noise/irregular component
    distvar.level=0*scale; % variance for level (mu)
    distvar.slope=2*scale; % variance for slope term (nu)
    
    distvar.acc=0*scale*(scale^4); % variance for acceleration term (xi)
    if strcmp(DataType,'GRACE-w-accelerations')
        distvar.slope=0;
        distvar.acc=2*scale;
    end
    distvar.cycle=[5 5]*scale; % variance for each cycle term (c)
    distvar.covcyclestar=distvar.cycle*0;
    
    distvar.AR=1;
    % covariances
    distvar.covirr=distvar.irr*.1;
    distvar.covslope=distvar.slope*0.1;
    
elseif strcmp(DataType,'TideGauge')
    scale=1;
    distvar.level=0.0;
    distvar.irr=300;
    distvar.acc=0;
    distvar.slope=0.0034;
    distvar.cycle=[1,0.5;1,0.5]; % variance for each cycle term
    distvar.AR = [1000]; % Variance for AR processes. Is optimized by EM algoithm if settings.fixprocvar=false
elseif strcmp(DataType,'Nile')
    scale=1; % was scale = 10
    %distvar.irr=15099;
    distvar.irr=10^4.2;
    %distvar.level=1469.1*scale;
    distvar.level=10^3.1;
    distvar.acc=0;
    distvar.slope=[];
    distvar.cycle=[];
    distvar.covcycle=[];
elseif strcmp(DataType,'InSAR')
    scale=10;
    distvar.irr=1;
    distvar.level=0;
    distvar.acc=0;
    distvar.slope=1;
    distvar.cycle=[1]; % variance for each cycle term
elseif strcmp(DataType,'InSARHarlingen1') || strcmp(DataType,'InSARHarlingen2')
    distvar.irr=1;
    distvar.level=0;
    distvar.acc=0;
    distvar.slope=1;
    distvar.cycle=[1]; % variance for each cycle term
elseif strcmp(DataType,'synthetic')
    scale=1;
    distvar.irr=100;
    distvar.level=0*scale;
    distvar.acc=0;
    distvar.slope=1;
    distvar.cycle=[1;1]; % variance for each cycle term
    % set parameters
    syn.cycle = [50 0];% c / c* recursive formulation
    syn.mu(1)=0;
    syn.nu(1)=0.;
    syn.xi(1)=0;
    synvar.level=0;
    synvar.slope=0.01;
    synvar.acc=0;
    synvar.cycle=[1 0];
    synvar.irr=10;
elseif strcmp(DataType,'KSI_UK')
    scale=1;
    distvar.irr=3e-3;
    distvar.level=9e-4;
    distvar.cycle=5e-7;
    distvar.covcyclestar=distvar.cycle*0;
elseif strcmp(DataType,'Paper_Willen_etal_2021_basin_21')
    distvar.irr=0.1;
    distvar.level=0; % variance for level (mu)
    distvar.slope=1e-1; % variance for slope (nu)
    distvar.acc=0; % variance for acceleration term (xi)
    distvar.cycle=[.1 .1]; % cycle variances, one for each period
   
    distvar.covcyclestar=distvar.cycle*0;% auxiliary cycle term covariances
    distvar.AR=[1e0]; % variance autoregressive component
    
    % covariances
    
    distvar.covAR=distvar.AR*0;% no correlations between AR terms
    % use ice density to relate variances of ALT-FDM and GRACE-cSMBA
    runsettings.icedensity = 0.917;
    
    runsettings.scalevariance= 1/runsettings.icedensity^2;
    
    distvar.slope(2)=distvar.slope(1)/runsettings.scalevariance; % GRACE-cSMBA has slightly higher variance
    
    distvar.covslope=prod(sqrt(distvar.slope))*0;% fully correlated
    
end

clear scale


%% Load time series

if ~strcmp(DataType,'synthetic')
    Temp=importdata(DataFile);
    if settings.regression
        Temp2=importdata(RegressorFile);
    end
    
    if strcmp(DataType,'GRACE') || strcmp(DataType,'GRACE-w-accelerations')
        year=Temp(:,1);
        month=Temp(:,2);
        y=Temp(:,3);
        % time
        time=year+(month-0.5)/12;
        ntimes=length(time);
        clear year month
    elseif strcmp(DataType,'TideGauge')
        time=Temp.meas_time/365.24;
        y=Temp.meas_height;
        ntimes=length(time);
        if settings.regression
            tseries.regressors(1,:)=Temp2.var_uwind;
            tseries.regressors(2,:)=Temp2.var_pressure;
        end
    elseif strcmp(DataType,'Nile')
        time=Temp(:,1);
        ntimes=length(time);
        y=Temp(:,2);
    elseif strcmp(DataType,'InSAR')
        time=Temp.t;
        ntimes=length(time);
        y=Temp.vert';
    elseif strcmp(DataType,'InSARHarlingen1') || strcmp(DataType,'InSARHarlingen2')
        time=Temp.time';
        ntimes=length(time);
        y=Temp.serie';
    elseif strcmp(DataType,'KSI_UK')
        time=Temp.time;
        ntimes = length(time);
        y=log(table2array(Temp.data));
        % remove bias
        y=y-mean(y,'omitnan');
    elseif strcmp(DataType,'Paper_Willen_etal_2021_basin_21')
        time=Temp.data(:,1);
        ntimes=length(time);
        y(1,:)=Temp.data(:,17)-Temp.data(:,18);
        y(2,:)=Temp.data(:,4);
       tseries.name{1}='altimetry-FDM';
       tseries.name{2}='GRACE-SMBA';
       tseries.generalname='basin 21';
       
       figure;plot(time,y(1,:))
       yyaxis right
       plot(time,y(2,:))
       legend(tseries.name{1},tseries.name{2})
    end
    clear Temp Temp2
else
    %% create synthetic data set
    
    [time,t,y,ntimes,syn]=synthetictimeseries(syn,synvar,settings);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDITING TEXT BELOW THIS LINE AT YOUR OWN RISK  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% put all time series information in structure
%    tseries.name = strcat(settings.batchname,' ',num2str(nsample));
if ~isfield(tseries,'name')
tseries.name = strcat(settings.batchname);
end
tseries.time=time;
tseries.t=[];
tseries.y=y;
tseries.Y=[];
tseries.ntimes=ntimes;

clear y time ntimes

%% Deterministic fitting using least squares
if settings.lsq
    % estimate using least squares
    [estlsq]=LeastSquaresTrend(tseries,settings);
    
    
    estlsq.residualRMS;
    % make plots
    figure;subplot(2,1,1)
    hold on;
    plot(tseries.time,tseries.y,'k',tseries.time,estlsq.fit)
    plot(tseries.time,estlsq.trend)
    legendstring={'observation','fit'};
    if settings.acclsq
        legendstring=[legendstring,{'quadratic trend'}];
    else
        legendstring=[legendstring,{'linear trend'}];
    end
    
    xlabel('time');ylabel('observation and fit')
    title(strcat({'least squares estimate for '},{tseries.name}))
    legend(legendstring);
    
    subplot(2,1,2); hold on
    plot(tseries.time,tseries.y,'.-k')
    plot(tseries.time,estlsq.residuals)
    xlabel('time');ylabel('residuals')
    legend('observation','residuals')
    
    %% Spectral analysis of original data and estimates using least squares
    % (Average) time between measurements
    Period=(tseries.time(end)-tseries.time(1))/(tseries.ntimes-1);
    if settings.continuoustime
        % do Lomb-Scargle
        [ypsd,f]=fastlomb(tseries.y,tseries.time);
        [y_fitLQ_psd,f_fit]=fastlomb(estlsq.fit,tseries.time);
        [TrendLQ_psd,f_trend]=fastlomb(estlsq.trend,tseries.time);
        [ResidLQ_psd,f_resid]=fastlomb(estlsq.residuals,tseries.time);
        if settings.cycle
            for i=1:length(settings.periodscycle)
                [CycleLQ_psd{i},f_Cycle{i}]=fastlomb(estlsq.cycles(i,:),tseries.time);
            end
        end
    else
        
        % Sampling frequency (average in case of inequally space time)
        Frequency=1/Period;
        % do Fourier analysis
        [f,ypsd]=SpectralAnalysis(tseries.y,tseries.time,Frequency);
        [f_fit,y_fitLQ_psd]=SpectralAnalysis(estlsq.fit,tseries.time,Frequency);
        [f_trend,TrendLQ_psd]=SpectralAnalysis(estlsq.trend,tseries.time,Frequency);
        [f_resid,ResidLQ_psd]=SpectralAnalysis(estlsq.residuals,tseries.time,Frequency);
        if settings.cycle
            for i=1:length(settings.periodscycle)
                [f_Cycle{i},CycleLQ_psd{i}]=SpectralAnalysis(estlsq.cycles(i,:),tseries.time,Frequency);
            end
        end
    end
    
    % Plot single-sided amplitude spectrum.
    figure
    legendstring={'observation','fit','trend','residual'};
    semilogx(f,ypsd,'k') ; hold on
    semilogx(f_fit,y_fitLQ_psd)
    semilogx(f_trend,TrendLQ_psd)
    semilogx(f_resid,ResidLQ_psd)
    if settings.cycle
        for i=1:length(settings.periodscycle)
            semilogx(f_Cycle{i},CycleLQ_psd{i});
            legendstring=[legendstring,{strcat('cycle',num2str(i))}];
        end
    end
    title(strcat({'Single-Sided Amplitude Spectrum of y(t) using least squares '},{tseries.name}))
    
    xlabel('Frequency (1/yr)')
    ylabel('|Y(f)|')
    xlim([min(1/(tseries.ntimes*Period)) max(f)])
    legend(legendstring);
end


%% run state space analysis %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tseries,model,distvar,settings,estfilter,fitfilter,estsmooth,fitsmooth,varsmooth,ConvError,covsmooth,meanslope]=statespaceanalysis(tseries,settings,distvar);

%% plotting results %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot filtered estimate
% Plot one step ahead prediction (filtered estimate)
if settings.checkplots
    PlotFit(settings,tseries,fitfilter,[],model,estfilter,'forward filtered estimate ')
end

%% Plot smoothed estimate
PlotFit(settings,tseries,fitsmooth,varsmooth,model,estsmooth,'smoothed estimate ')


%% plot amplitude and phase information of cycle terms
PlotCycles(settings,tseries,fitsmooth,varsmooth)


%% plot autocorrelation state trend
PlotCovariance(tseries,model,covsmooth,estsmooth,settings,'autocorrelation')

%% determine and print averaged slope from time variable slope
PlotMeanSlope(tseries,settings,fitsmooth,meanslope,varsmooth,model)

%% Spectral analysis of original data and state space estimates
[psd]=SpectralPlot(settings,tseries,estsmooth,fitsmooth);


%% statistical test of irregular term for normality
%Single sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
% this is to test the gaussian distribution of the irregular term
for p=1:tseries.ntseries
    residual=estsmooth.epsilon(p,:);
    H=kstest((residual-mean(residual))/std(residual));
    if H==0
        disp('do not reject the null hypothesis of gaussian residuals of the residual at 5% significance')
    else
        disp('Reject the null hypothesis of gaussian residuals at the 5% significance level.')
    end
    
    lag=round(length(tseries.time)-1);
    significancelevel=0.05;
    figure
    [H, pValue, Qstat, CriticalValue]=lbqtest1(residual,lag,significancelevel);
    title(strcat('Residual Autocorrelation function, time serie nr:',num2str(p)))
end
%% Ljung-Box Q-statistic lack-of-fit hypothesis test of AR
if settings.AR
    for p=1:tseries.ntseries
        figure
        lbqtest1(fitsmooth.AR{p},lag,significancelevel)
        title(strcat('AR(1) Autocorrelation function',num2str(p)))
    end
end

%% Diagnostic test
% Standardized innovations
for p=1:tseries.ntseries
    if tseries.ntseries==1
      v=estfilter.v;
      F=estfilter.F;
    else
    v=squeeze(estfilter.v(p,:))';
    F=sqrt(squeeze(estfilter.F(p,p,:)));
    end
norm_res=v./F; %
%% QQ-Plot of standardized innovations
% ====================================
figure; qqplot(norm_res)
title('QQ Plot of Standardized Innovations');

%% Test the standardised innovations for outliers
%%=========================================================================

idx = ~isnan(norm_res);
norm_res = norm_res(idx); N=length(norm_res); %exclude nan data values
figure;
plot(tseries.time(idx),norm_res,'ro');hold on;
xlabel('time [years]');ylim([-3,3]);
title('Standardised innovations and 95% confidence interval');
plot(tseries.time(idx),ones(1,N)*1.96,'-k');
plot(tseries.time(idx),-ones(1,N)*1.96,'-k');
legend('Standardized innovations','95% confidence interval');
hold off;

%% Test on independence of innovations using the Box-Ljung statistic
% ==================================================================
maxLag=round( length(tseries.time(idx))/10 );
maxLag=min(maxLag,20);
significancelevel=0.05; DoF = maxLag - model.dim+1;
[H,pValue,Qstat,CriticalValue]=lbqtest1(norm_res,maxLag,significancelevel,DoF);

fprintf('\nBox-Ljung test on independence of innovations\n')
if H==0
    disp('Do not reject the null hypothesis')
else
    disp('Reject the null hypothesis')
end
fprintf('Significance level is %4.2f\n',significancelevel);
fprintf('p value is %4.2f; teststatistic is %4.2f; critical value is %4.2f\n\n',pValue,Qstat,CriticalValue);
title('Autocorrelation sequence of innovations');xlabel('lag');


%% Test on homoscedasticity of innovations using Goldberg-Quandt test
% =====================================================================
siglevel=0.05; DoF = maxLag - model.dim+1;
h = round(N/3);
Cvalue = finv(1-siglevel,h,h);
rr = norm_res - mean(norm_res);
stat = sum(rr(N-h+1:end).^2)/sum(rr(1:h).^2);
if stat < 1
    stat = 1/stat;
end
fprintf('Goldfeld-Quandt test on homoscedasticity of innovations\n');
if stat < Cvalue
    fprintf('Do not reject the null-hypothesis; test statistic = %8.3f\n',stat);
else
    fprintf('Reject the null hypothesis; test statistic = %8.3f\n',stat);
end
fprintf('One-sided F-test; significance level = %8.2f; critical value = %8.5f\n\n',siglevel,Cvalue);



%% Bowman Shenton test on normality of innovations
% ===============================================
siglevel = 0.05; Cvalue = chi2inv(1-siglevel,2);
rr = norm_res - mean(norm_res);
S = skewness(rr); K = kurtosis(rr);
stat = N*(S^2/6 + (K-3)^2/24);
fprintf('Bowman-Shenton test on normality of innovations\n');
if stat < Cvalue
    disp('Do not reject the null hypothesis')
else
    disp('Reject the null hypothesis')
end
fprintf('Significance level is %4.2f\n',siglevel);
fprintf('teststatistic is %4.2f; critical value is %4.2f\n',stat,Cvalue);


%% Test the smoothed residuals on normality using Anderson-Darling test
% ====================================================================
residual=estsmooth.epsilon(p,:);
[H,pp,adstat,cv]=adtest(residual,'Alpha',siglevel);
fprintf('\nAnderson-Darling test on normality of Kalman smoother residuals\n');
if H==0
    disp('Do not reject the null hypothesis');
else
    disp('Reject the null hypothesis');
end
fprintf('Significance level is %4.2f\n',siglevel);
fprintf('p value is %4.2f; teststatistic is %4.2f; critical value is %4.2f\n',pp,adstat,cv);



end

%% plot difference between kalman & LS
if settings.lsq
    for p=1:tseries.ntseries
        figure;hold on
        subplot(2,1,1);hold on;clear handle
        handle(1)=plot(tseries.time,tseries.Y,'k');
        handle(2)=plot(tseries.time,fitsmooth.fit{p});
        handle(3)=plot(tseries.time,estlsq.fit);
        
        handle(4)=plot(tseries.time,fitsmooth.trend{p});
        handle(5)=plot(tseries.time,estlsq.trend);
        legend(handle,'data','state space fit','least squares fit','state space trend','least squares trend')
        
        title('least squares fit vs. kalman fit')
        subplot(2,1,2);hold on;plot(tseries.time,abs(estlsq.residuals),tseries.time,abs(estsmooth.epsilon))
        legend('least squares residuals','state space residuals')
        xlabel('time')
        title('absolute residuals')
    end
end



%% in case of synthetic run
if strcmp(DataType,'synthetic')
    
    figure; hold on
    subplot(2,1,1);hold on
    plot(tseries.time,tseries.y,'k')
    plot(tseries.time,fitsmooth.fit{1},'r');
    plot(tseries.time,estsmooth.epsilon,'c')
    title('original time series and fit')
    legend('time series','fit','residual')
    
    subplot(2,1,2); hold on
    plot(tseries.time,syn.mu(1:tseries.ntimes),'g')
    
    plot(tseries.time,syn.cycle(1:tseries.ntimes),'b')
    plot(tseries.time,fitsmooth.trend{1},'go')
    if settings.cycle
        for i=1:settings.numbercycles
            plot(tseries.time,fitsmooth.cycle{1}(i,:),'bo')
        end
    end
    title('original components and fit')
    legend('trend','cycles')
end

%% save results
if settings.saveresults
    if isfield(tseries,'name')
        savefile=strcat(settings.maindir,'/',settings.savedir,'/savedresults_',tseries.generalname,'.mat');
    else
        savefile=strcat(settings.maindir,'/',settings.savedir,'/savedresults.mat');
    end
    
    save(savefile,'settings','tseries','fitsmooth','psd','varsmooth','estsmooth','covsmooth','ConvError','distvar','estfilter','meanslope')
    
    % save as ascii
    WriteStateSpaceAscii(settings,tseries,model,fitsmooth,estsmooth,varsmooth,distvar)
    
    
end

%end


