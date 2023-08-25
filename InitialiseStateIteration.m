function [estfilter] = InitialiseStateIteration(tseries,settings,model,distvar)
% InitialiseState.m provides initial values for the filtered state and
% variance
% 
%
% HOW:     [estfilter] = InitialiseStateIteration(tseries,settings,model)
%
% Input:    
%       tseries, structure of which this function uses:
%           tseries.t [epochs x 1]          normalised time vector
%           tseries.Y                       observation vector
%           tseries.ntimes                  number of epochs
%           tseries.regressors [NRegr x epochs] time series functioning as explanatory variable
%           tseries.missingepochs           indexes of 'time' with missing data
%
%       settings, structure of which this function uses:
%           settings.slope                  determines whether a slope will be estimated
%           settings.acc                    determines whether accelerations will be estimated
%           settings.cycle                  determines whether cycle components will be estimated
%           settings.regression             determines whether regression is used (using user supplied regressors)
%           settings.intervention           determines whether step intervention effects will be estimated
%           settings.numbercycles           number of cycle terms
%           settings.lambda                 2pi/period of cycles
%           settings.continuoustime         whether time is non-integer (thus continuous)
%           settings.type                   type of formulation of the state space equations
%           settings.tau [N_interventions]  index of time of the intervention (index of normalized time)
%           settings.filter                 type of Kalman formulation
%           settings.maxiterinit            number of iterations for initialisation of a0 and P0 filter estimates
%           settings.largeP                 large number for initialisation of filtered state error variance P
%                               
%       model structure
%           model.dim (n)               dimension of the state vector
%
% Output: 
%       estfilter, structure of which this function outputs
%           estfilter.a [n x epochs]         initial filtered state
%           estfilter.P [n x n  x epochs]    initial state variance
% 
%
% Note:    
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
% uses: Kalman,Smoother
% 
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 june 2014 DBT Broerse
%
% Version 1.1 August 2014 DBT Broerse
% uses a Kalman filter-smoother iteration to set a(1) and P(1)
%
% Version 1.2 September 2014 DBT Broerse
% - added regression
% - added step intervention
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
%
%----------------------------------------------------------------------------
% remarks: still to be added: diffuse initialisation
%----------------------------------------------------------------------------

%
% start program

% initialisation and assume zero state prediction and large variance for t=1
% diffuse initialisation may be considered later on




%% initialise for first time with very high error variance
% prediction for state vector
estfilter.a=zeros(model.dim,tseries.ntimes);
% state variance
estfilter.P=zeros(model.dim,model.dim,tseries.ntimes);
estfilter.P(:,:,1)=eye(model.dim)*settings.largeP;

% exceptions for when initial state is known


% set AR to zero at first epoch
if settings.AR
if settings.startARat0 
    for i=1:length(model.indexAR)
    estfilter.P(model.indexAR{i},model.indexAR{i},1)=0;
     %   estfilter.P(model.indexAR{i},model.indexAR{i},1)=distvar.AR(i);
    end
end
end

% add Q to P, this has no effect for those states with unknown initial
% state, but adds variances and covariance for states with known initial
% value
%disp('adding Q1 to P1')
estfilter.P(:,:,1)=estfilter.P(:,:,1)+model.R*model.Q(:,:,1)*model.R';





end