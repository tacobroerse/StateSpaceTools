function [settings,distvar]=SetDefaults
% set defaults for state space analysis
%
% HOW:     [settings,distvar]=SetDefaults
%
% Output:
%           settings
%           distvar             disturbance variances of process and irregular component
%----------------------------------------------------------------------------
% uses: none
% 
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 October 2014 DBT Broerse
%
% Version 1.1 September 2017 DBT Broerse
% added settings for saving plots and results
%
% Version 1.2 June 2018 DBT Broerse
% added AR settings
%
%----------------------------------------------------------------------------
% remarks: 
%----------------------------------------------------------------------------

% general settings
settings.type = 'Normal_Slope_and_Acc'; % type of state space formulation ('Normal_Slope_and_Acc' is default, you may use StateSpaceType=[])
settings.filter = 'DurbinKoopman'; % kalman filter and smoother type (default is 'DurbinKoopman', other option SquareRoot)
settings.smoother = []; % kalman filter and smoother type (default is 'DurbinKoopman', other option SquareRoot)
settings.distsmoother = [];
settings.multivariate = 0; % 0 : univariate ; 1 : multivariate
settings.optvar = 'EM'; % variance optimizer
settings.checktime = false; % check whether time is integer
settings.continuoustime = true; % if false time steps are equal
settings.fixprocvar = false; % if true no variances for the irregular component and process noise are estimated
settings.checkplots = false; % produce extra plots for checks 
settings.lsq = false; % do a least squares analysis for comparison
settings.maxiterEM= 1000; % Iteration maximum of EM algorithm
settings.convEM = 1e-3; % Convergence criterium EM algorithm (min difference in log likelihood between two consequtive iterations)
settings.plotEM = true; % making plots of convergence of variance parameters in EM algorithm
settings.MaxRelDeviationLogL=0;% maximum relative allowed deviation in EM for log likelihood
settings.maxiterinit=2; % number of iterations for initialisation of a0 and P0 filter estimates
settings.largeP=1e5; % large number for initialisation of filtered state error variance P
settings.normalise = true; % normalise time
settings.normalisetimeby='median';% method of choosing characteristic frequency : mean or median
settings.estimatemeanslope = false; % estimate mean slope

settings.savedir='results/'; % directory for saving
settings.saveplots=1;% whether or not to save resulting plots
settings.saveresults=0;

settings.plotconstrainperiod=0;% whether or not to show periods without data, at the beginning and end of timeseries


%% Set which components should be included in the state space model (and thus estimated)
settings.slope = true; %slope, for a local linear trend model a.k.a. integrated random walk if stochastic (process variance>0), linear drift if deterministic (process variance=0), keep at true for smooth trends
settings.acc = false; %acceleration, however seems not so useful when trend is stochastic
settings.acclsq = false; % acceleration for least squares
settings.cycle = false; %cycle terms (harmonic series)
settings.regression = false; % regression (fitting functions)
settings.intervention = false; % step function at time tau
settings.numbercycles = 0;  %number of cycle terms
settings.numberinterventions = 0;   %number of interventions
settings.numberregressors = 0; %number of regressors   
settings.lambda  = [];  %normalised cycle frequency
settings.tau  = [];   %time index, startpoint for step function
settings.periodscycle = []; 
settings.interventiontime=[];
settings.AR = false;
settings.ARphi   = []; % Set autoregressive coefficients
settings.optAR   = false; % Optimize choice of autoregressive component. Use only for order-1 AR processes
settings.ARorder = 1;  % AR order, only 1 is currently implemented
settings.ARphiRange=[0:0.1:1];  % range for searching optimal phi value AR process
settings.ARphiMulti = 0; % separate AR phi value for each time series
settings.startARat0= 0;% require AR process to start at 0 at first time step
% time series structure
%tseries.regressors=[];% time series of regressors
%tseries.missingepochs=[];% index of missing epochs (only for integer time = equal time steps)
settings.removeinterventionepochs = 0; % remove epochs coinciding with interventions

% initial variances
distvar.irr=1; % variance for noise/irregular component
distvar.level=0; % variance for level (mu)
distvar.slope=1; % variance for slope term (nu)
distvar.acc=0; % variance for acceleration term (xi)
distvar.cycle=[]; % variance for each cycle term (c)
distvar.AR=[];
distvar.covcyclestar=[];
settings.distvarcovariancecycles=0;% whether or not covariance in cycle disturbance should be included (correlation = 1)
% following applies to multivariate models, whether or not to include
% covariance between different time series for separate components
settings.multivariatedistcovarianceirr=0;
settings.multivariatedistcovarianceslope=0;
settings.multivariatedistcovariancelevel=0;
settings.multivariatedistcovarianceacc=0;
settings.multivariatedistcovarianceAR=0;
%timevariableregressor=0;% if set to zero, a fixed regression coefficient will be estimated
end