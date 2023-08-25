function [fit]=PostProcessState(esttemp,settings,model,tseries)
% PostProcessState.m computes estimates of the several components of the time series
% based on filtered/smoothed estimates of the state vector.
% uses the same state space formulation as StateSpace.m. In case state space
% model is changed, PostProcessState.m should be adapted accordingly.
%
% HOW:     [fit]=PostProcessState(estfilter,settings,model,tseries)
%
% Input:
%       esttemp, structure that should contain
%           esttemp.alpha or esttemp.a      state vector
%           esttemp.epsilon or esttemp.v    irregurlar component or prediction error
%
%       tseries, structure of which this function uses
%           tseries.t                       normalised time vector
%           tseries.period                  time series sampling periodicity
%
%       settings, structure of which this function uses
%           settings.slope                  determines whether a slope will be estimated
%           settings.acc                    determines whether accelerations will be estimated
%           settings.cycle                  determines whether cycle components will be estimated
%           settings.regression             determines whether regression is used (using user supplied regressors)
%           settings.intervention           determines whether step intervention effects will be estimated
%           settings.numbercycles           number of cycle terms
%           settings.numberregressors       number of regressors
%           settings.numberinterventions    number of interventions
%           settings.lambda                 2pi/period of cycles
%           settings.AR                     determines whether AR components will be estimated

%           settings.type                   type of formulation of the state space equations
%                                           - 'Normal_Slope_and_Acc' default
%       model, structure of which this function uses
%           model.Z                         design matrix
%           model.index                       dimension of state vector
%           model.index[component]          indices of the several components
%
%
% Output:
%       fit, structure
%           fit.fit [epochs x 1]            estimated time series
%           fit.trend [epochs x 1]          trend component
%           fit.bias [1]                    bias (deterministic) or start point (stochastic)
%           fit.slope [epochs x 1]          time variable slope (unnormalized)
%           fit.acc [epochs x 1]            time variable acceleration (unnormalized)
%           fit.cycle {ncycl}[epochs x 1]   cycle components
%           fit.cycleAmpl [idem]            cycle component amplitudes
%           fit.cyclePhase [idem]           cycle component phases
%           fit.regression [n_reg x epochs] regression
%           fit.regressionCoef [n_reg]      regression coefficients (1 per regressor)
%           fit.intervention [n_pulse x ep] pulse intervention
%           fit.delta [n_pulse_intervent.]  change in level at time t(tau)
%           fit.distRMS                     RMS of disturbance component
%           fit.AR                          AR component
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
% uses: none
%
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 June 2014 DBT Broerse
%
% Version 1.1 August 2014 DBT Broerse
% - added state error variance and corresponding output
%
% Version 1.1 August 2014 DBT Broerse
% - added state error variance and corresponding output
% - denormalizing slope and acceleration
%
% Version 1.2 September 2014 DBT Broerse
% - added regression
% - added step intervention
% - removed state space options using lag operator
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
% - continuous time and integer time functions merged
%
% Version 2.1 March 2015 DBT Broerse
% - changed dimension of Z
% - add option for unnormalised slope and acceleration
%
% Version 2.2 June 2018 T Frederikse / DBT Broerse
% - added AR
%
% Version 2.3 June 2018 DBT Broerse
% - extended to multivariate models
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
%----------------------------------------------------------------------------
%error(nargchk(4,4,nargin));
narginchk(4,4);

% because this function should be able to handle both input from filtered
% estimates as well as smoothed estimates, special care is needed for the
% elements of structure est.

if any(strcmp('alpha',fieldnames(esttemp)))
    est.alpha=esttemp.alpha;
elseif any(strcmp('a',fieldnames(esttemp)))
    est.alpha=esttemp.a;
end

if any(strcmp('epsilon',fieldnames(esttemp)))
    est.epsilon=esttemp.epsilon;
elseif any(strcmp('v',fieldnames(esttemp)))
    est.epsilon=esttemp.v;
end

clear esttemp;

if isempty(settings.type);settings.type='Normal_Slope_and_Acc';end

if ~strcmp(settings.type,'Normal_Slope_and_Acc')
    error(strcat('state space type ',settings.type,' not supported, exiting'))
end


fit.fit = [];
fit.trend = [];
fit.bias = [];
fit.slope = [];
fit.acc = [];
fit.cycle = [];
fit.cycleampl = [];
fit.cyclephase = [];
fit.regressioncoef = [];
fit.delta = [];
fit.regression = [];
fit.intervention = [];
fit.AR = [];


% check whether Z is time dependent or not


switch settings.continuoustime
    
    case 0
        
        p=1;% no multivariate allowed, but added for consistency with continuous time
        % make mask
        MaskTrend = zeros(1,model.dim);
        MaskCycle = [];
        
        % there are a couple of issue that are not updated yet for integer
        % time
        p=1;
        % no time dependency in Z
        
       
        MaskTrend(model.indexlevel{p})=1;
        
        
        if settings.cycle
            for i=1:settings.numbercycles
               
                
                MaskCycle{i}=zeros(1,model.dim);
                MaskCycle{i}(model.indexcycles{p}((i-1)*2+1))=1;
                
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % create estimated time series
        
        fit.fit{p} = model.Z*est.alpha;
        
        %% create trend component
        
        fit.trend{p} = MaskTrend*est.alpha;
        
        %% create cycle components
        if settings.cycle
            for i=1:settings.numbercycles
                
       
                fit.cycle{p}(i,:)=MaskCycle{i}*est.alpha;
                % phase is modulus of atan(c_t*,c_t)-time*cycle_frequency
                j=(i-1)*2+1;
                % phase is modulus of atan(c_t*,c_t)-time*cycle_frequency
                
                fit.cyclephase{p}(i,:)=mod(-atan2(est.alpha(model.indexcycles{1}(j+1),:),est.alpha(model.indexcycles{1}(j),:))-tseries.t'*settings.lambda(i),2*pi);
                % amplitude is sqrt(c_t^2 + c_t*^2)
                
                fit.cycleampl{p}(i,:)=sqrt(est.alpha(model.indexcycles{1}(j+1),:).^2+est.alpha(model.indexcycles{1}(j),:).^2);
%                
%                 fit.cycle{i}=MaskCycle{i}*est.alpha;
%                 % phase is modulus of atan(c_t*,c_t)-time*cycle_frequency
%                 fit.cyclephase{i}=mod(-atan2(est.alpha(CycleStarIndex(i),:),est.alpha(CycleIndex(i),:))-(tseries.t)'*settings.lambda(i),2*pi);
%                 % amplitude is sqrt(c_t^2 + c_t*^2)
%                 fit.cycleampl{i}=sqrt(est.alpha(CycleStarIndex(i),:).^2+est.alpha(CycleIndex(i),:).^2);
            end
        end
        
        
        
        %% bias (or start position), slope and trend
        %
        % for now only possible for this state space formulation (since
        % acceleration and slope are not directly estimated in the formulation
        % using lag operators, may be done using numerical differentiation?
        
        if strcmp(settings.type,'Normal_Slope_and_Acc')
            % create bias
            fit.bias=est.alpha(model.indexlevel{p},1);
            if settings.slope
                % create slope
                fit.slope{p} = est.alpha(model.indexslope{p},:);
                % de normalise
                fit.slope{p} = fit.slope{p}/tseries.period;
                % create acceleration
                if settings.acc
                    fit.acc{p} = est.alpha(model.indexacc{p},:);
                    % de normalise
                    fit.acc{p} = fit.acc/tseries.period^2;
                end
            end
        end
        
        %% autoregressive
        if settings.AR
            fit.AR{p} = est.alpha(model.indexAR{p},:);
        end
        
        % continuous time
    case 1
        
        % time dependency in Z
        for p=1:tseries.ntseries
            % create masks
            
            MaskTrend{p} = zeros(1,model.dim);
            MaskTrend{p}(model.indexlevel{p})=1;
            
            if settings.cycle
                for i=1:settings.numbercycles
                    MaskCycle{p,i}=zeros(1,model.dim);
                    MaskCycle{p,i}(model.indexcycles{p}((i-1)*2+1))=1;
                end
            end
            
            if settings.regression
                for i=1:settings.numberregressors
                    MaskRegressor{p,i}(1:tseries.ntimes)=model.Z(p,model.indexregressors{p}(i),:);
                end
            end
            
            
            if settings.intervention
                for i=1:settings.numberinterventions
                    MaskIntervention{p,i}(1:tseries.ntimes)=model.Z(p,model.indexinterventions{p}(i),:);
                    
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % create estimated time series
            
            for i=1:tseries.ntimes
                fit.fit{p}(i) = model.Z(p,:,i)*est.alpha(:,i);
            end
            
            %% create trend component
            
            fit.trend{p} = MaskTrend{p}*est.alpha;
            
            %% create cycle components
            if settings.cycle
                for i=1:settings.numbercycles
                    
                    fit.cycle{p}(i,:)=MaskCycle{p,i}*est.alpha;
                    % phase is modulus of atan(c_t*,c_t)-time*cycle_frequency
                    j=(i-1)*2+1;
                    fit.cyclephase{p}(i,:)=mod(-atan2(est.alpha(model.indexcycles{p}(j+1),:),est.alpha(model.indexcycles{p}(j),:))-(tseries.t)'*settings.lambda(i),2*pi);
                    
                            disp(' please review definition phase')
                         
                    % amplitude is sqrt(c_t^2 + c_t*^2)
                    
                    fit.cycleampl{p}(i,:)=sqrt(est.alpha(model.indexcycles{p}(j+1),:).^2+est.alpha(model.indexcycles{p}(j),:).^2);
                    
                   
                end
            end
            
            %% bias (or start position), slope and trend
            
            % create bias
            fit.bias{p}=est.alpha(model.indexlevel{p},1);
            if settings.slope
                % create slope
                fit.slope{p} = est.alpha(model.indexslope{p},:);
                if settings.normalise
                    % de normalise
                    fit.slope{p} = fit.slope{p}/tseries.period;
                    %disp('denormalise slope')
                end
                % create acceleration
                if settings.acc
                    fit.acc{p} = est.alpha(model.indexacc{p},:);
                    % de normalise
                    if settings.normalise
                        fit.acc{p} = fit.acc{p}/tseries.period^2;
                        %    disp('denormalise acceleration')
                    end
                end
            end
            
            if settings.regression
                for i=1:settings.numberregressors
                    fit.regressioncoef{p}(i)=est.alpha(model.indexregressors{p}(i));
                    fit.regression{p}(i,:)=MaskRegressor{p,i}*fit.regressioncoef{p}(i);
                end
            end
            
            
            if settings.intervention
                for i=1:settings.numberinterventions
                    fit.delta{p}(i)=est.alpha(model.indexinterventions{p}(i));
                    fit.intervention{p}(i,:)=MaskIntervention{p,i}*fit.delta{p}(i);
                end
            end
            
            if settings.AR
                fit.AR{p} = est.alpha(model.indexAR{p},:);
            end
        end
        
        
end

for p=1:tseries.ntseries
    % rms of irregular
    
    fit.irregularRMS{p}=sqrt(1/length(tseries.t)*sum(est.epsilon(p,:).^2));
    if settings.AR
        fit.ARRMS{p}=sqrt(1/length(tseries.t)*sum(fit.AR{p}.^2));
        % total irregular and AR
        fit.residualRMS=sqrt(fit.ARRMS{p}.^2+fit.irregularRMS{p}^2);
    end
end

end

