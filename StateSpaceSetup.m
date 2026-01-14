function [model,distvar] = StateSpaceSetup(tseries,settings,distvar)
% StateSpace.m determines matrices needed by the Kalman filter
%
% HOW:     [model,distvar] = StateSpaceSetup(tseries,settings,distvar)
%
% Input:
%       tseries, structure of which this function uses:
%           tseries.t [epochs x 1]          normalised time vector
%           tseries.ntimes                  number of epochs
%           tseries.regressors [NRegr x epochs] time series functioning as explanatory variable
%           tseries.ntseries                number of time series
%
%       settings, structure of which this function uses:
%           settings.slope                  determines whether a slope will be estimated
%           settings.acc                    determines whether accelerations will be estimated
%           settings.cycle                  determines whether cycle components will be estimated
%           settings.regression             determines whether regression is used (using user supplied regressors)
%           settings.intervention           determines whether step intervention effects will be estimated
%           settings.numbercycles           number of cycle terms
%           settings.numberregressors       number of regressors
%           settings.numberinterventions    number of interventions
%           settings.lambda                 2pi/period of cycles
%           settings.AR                     determines whether an AR process will be estimated
%           settings.continuoustime         whether time is non-integer (thus continuous)
%           settings.type                   type of formulation of the state space equations
%                                           - 'Normal_Slope_and_Acc' default
%           settings.tau [N_interventions]  index of time of the intervention (index of normalized time)
%
%
%       distvar, structure with:
%           distvar.level                       variance of level term (usually zero to create smooth trend)
%           distvar.slope                       variance of slope term (zeta in Durbin&Koopman)
%           distvar.acc                         variance of acceleration
%           distvar.cycle [settings.numbercycles x 1]    variance of cycle terms
%           distvar.covcyclestar                covariance between cycle and cycle star terms
%           distvar.irr                         variance of disturbance
%           distvar.AR                          variance of AR process
%
%           distvar.varcovirr [p x p]                  variance-covariance matrix irregular (for multivariate models)
%           distvar.covlevel  [p x p]                   covariance matrix level (for multivariate models)
%           distvar.covslope  [p x p]                   covariance matrix slope (for multivariate models)
%           distvar.covacc    [p x p]                  covariance matrix accelerations (for multivariate models)
%           distvar.covAR    [p x p]                  covariance matrix AR (for multivariate models)
%
%
% Output:   model, structure, which contains
%           model.Z   [p x n x epochs]        design matrix
%           model.T   [n x n x epochs]        transition matrix
%           model.Q   [n x n x epochs]        process noise matrix
%           model.R   [n x n]                 selection matrix
%           model.H   [p x p x epochs]                disturbance variance matrix
%           model.index                       indices of the state vector of the different components
%
%           distvar.eta                       disturbance variance vector
%
% where p is the number of time series (p = 1 for univariate models)
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
% note: preprocess and check distvar before using StateSpaceSetup using ProcessNoiseCheck
%
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 June 2014 DBT Broerse
%
% Version 1.1 September 2014 DBT Broerse
% - added regression
% - added step intervention
% - design matrix now time dependent
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
% - continuous time and integer time functions merged
%
% Version 2.1 March 2015 DBT Broerse
% - removed first dimension of Z matrix
% - added dt in definition of H matrix for continuous time
% - phased out state space definitions based on lag operators
%
% Version 2.2 April 2015 DBT Broerse
% - process noise vector eta is defined as well (for use in EM later on)
% - extra output
%
% Version 2.3 June 2018 DBT Broerse/T Frederikse
% - AR process included
%
% Version 2.4 June 2018 DBT Broerse
% - extension to multivariate models
% - made Z consistent with other matrices in terms of index order
%
% Version 2.5 June 2024 DBT Broerse
% - removed time step dependence from H for continuous time
%
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------
%
%%%% THEORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The basic time series model with time t=1,2..T is defined in state space form
% (Durbin & Koopman, section 3.2):
%
% y_t = mu_t + c_t + epsilon_t      epsilon ~ N(0,sigma_eps^2)      (eq.1)
%
% y_t is the observation
% mu_t is the slowly varying trend
% c_t is the cycle term (harmonic term)
% epsilon is irregular component called disturbance
%
% the observation equation is:
%
% y_t = Z_t * alpha_t + epsilon
% where alpha is the state vector, Z is the design matrix
%
% and the state equation is:
% alpha_t+1 = T_t * alpha_t + Q
% where T is the transition matrix, Q process noise matrix, often an
% identity matrix.
%
%%%% Trend part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mu_t+1 = mu_t + nu_t*dt + zeta_t  zeta_t  ~ N(0,sigma_zeta^2)
% nu_t+1 = nu_t + xi_t              xi_t  ~ N(0,sigma_xi^2)
%
% where dt = t(i+1) - t(i)
%
% Setting sigma_zeta^2=sigma_xi^2=0, we get a deterministic trend
% to make the trend more smooth, we set sigma_zeta^2 = 0. The resulting
% trend is called integrated random walk or local linear trend.
%
%
%
%%%% Quadratic part / acceleration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Including an acceleration term the trend becomes
%
%
% mu_t+1 = mu_t + nu_t*dt + 1/2*a_t*dt^2     zeta_t  zeta_t  ~ N(0,sigma_zeta^2 dt)
% nu_t+1 = nu_t + a_t*dt + xi_t              xi_t  ~ N(0,sigma_xi^2 dt)
% a_t+1  = a_t + chi_t                       chi_t ~ N(0,sigma_chi^2 dt)
%
%
% Setting sigma_zeta^2=sigma_xi^2=sigma_chi^2=0, we get a deterministic trend
% to make the trend more smooth, we set sigma_zeta^2 = 0. Including an
% acceleration term is not common and probably of little use in a
% stochastic setting since accelerations are already accomodated in the
% nu_t term.
%
%%%% Cycle part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% c_t = c cos(lambda t) + c* sin(lambda t)
%
% where we have a c_t for every cycle period
%
% using the addition rules for cosine and sine one can also write
% (Durbin and Koopman 3.2.4):
%
% c_t+1  = c_t cos(lambda*dt)  + c*_t sin(lambda*dt) + omega_t           (eq.4)
% c*_t+1 = -c_t sin(lambda*dt) + c*_t cos(lambda*dt) + omega_t 
%
% where omega_t ~ N(0,sigma_omega^2 dt)
%
% Continuous time version can be found in Harvey, equation 9.2.5, omitting
% the rho damping factor.
%
%%%% Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To find the scaling between y and regressor x_t
%
% Sum_j=1^k Beta * x_jt  for k regressors x_t
%
% with no process noise included
%
%%%% Step Intervention %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To find step function delta at time t(tau)
%
% delta * w_t
%
% with w_t = 0 for t < t(tau)
%      w_t = 1 for t => t(tau)
%
% with no process noise included
%
%%%% AR process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% see Frederikse et al 2016
%
% zeta_t+1 = zeta_t * phi^dt + psi 
% where psi ~ N(0,sigma_psi^2 dt)
%
% where dt = t(i+1) - t(i)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the state vector becomes:
%
% alpha = [mu_t nu_t] for linear trend
%
% or
%
% alpha = [mu_t nu_t a_t] for linear trend plus quadratic
%
%
% alpha = [mu_t nu_t a_t c1_t c1*_t c2_t c2*_t] for trend
% plus quadratic and two cycle terms
%
%
% alpha = [mu_t nu_t a_t c1_t c1*_t c2_t c2*_t Beta delta] for trend
% plus quadratic and two cycle terms, one regressor and one intervention
%
% etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% note: for integer time steps, i.e. settings.continuoustime is false
% dt is replaced by 1. which leads to time independent matrices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start program
dt=tseries.dt;
% different code for continuous time (general case) and integer time steps
switch settings.continuoustime
    
    case 1
        
        % set default
        if isempty(settings.type);settings.type='Normal_Slope_and_Acc';end
        % checks
        if settings.regression
            [nrow ncol]=size(tseries.regressors);
            if ncol~=tseries.ntimes
                disp('regressor has unequal length compared to observation vector')
                settings.ntimes
                ncol
                return
            end
        end
        
        % continuous time version
        if strcmp(settings.type,'Normal_Slope_and_Acc')
            %% Determine dimension of matrices
            model.indexslope=[];
            model.indexacc=[];
            model.indexcycles=[];
            model.indexregressors=[];
            model.indexinterventions=[];
            model.indexAR=[];
            % trend uses 2 dimensions (mu + v)
            % acceleration adds 1 extra
            % each cycle 2
            % each regressor 1
            % each intervention 1
            
            n=0; % current dimension count
            
            % loop on number of time series
            for p=1:tseries.ntseries
                n=n+1;
                % NOTE: for future changes, the dimensions of level, slope
                % and acceleration have to be sequential
                model.indexlevel{p} = n;
                
                if settings.slope
                    n=n+1;
                    model.indexslope{p} = n;
                    
                    if settings.acc
                        n=n+1;
                        model.indexacc{p} = n;
                    end
                end
                
                if settings.cycle
                    model.indexcycles{p}=[n+1:n+2*settings.numbercycles];
                    %   model.indexcycles=[model.indexcycles n+1:n+2*settings.numbercycles];
                    n=n+2*settings.numbercycles;
                end
                
                
                
                
                
                if settings.regression
                    model.indexregressors{p}=[n+1:n+settings.numberregressors];
                    %  model.indexregressors=[model.indexregressors n+1:n+settings.numberregressors];
                    n=n+settings.numberregressors;
                end
                
                if settings.intervention
                    model.indexinterventions{p}=[n+1:n+settings.numberinterventions];
                    %      model.indexinterventions=[model.indexinterventions n+1:n+settings.numberinterventions];
                    n=n+settings.numberinterventions;
                end
                
                if settings.AR
                    n=n+1;
                    model.indexAR{p}=n;
                    % add extra indices for higher order AR processes
                    if settings.ARorder > 1
                        n=n+settings.ARorder-1;
                    end
                end
                
                
            end
            
            
            
            %% make matrices
            
            % n : model dimension
            
           
                Z = zeros(tseries.ntseries,n,tseries.ntimes);
                T = zeros(n,n,tseries.ntimes);
                Q = zeros(n,n,tseries.ntimes);
           
            distvar.eta = zeros(1,n);
            
            %% H: disturbance variance
           
            if settings.multivariate
                
                for itimes=1:length(tseries.dt)
                    H(1:tseries.ntseries,1:tseries.ntseries,itimes) = distvar.varcovirr;
                end
            else
                %H = diag(dt)*distvar.irr; % dt should not enter equation
                H = eye(length(dt))*distvar.irr;
            end
            
            
            %% Z: design matrix
            
            for p=1:tseries.ntseries
                Z(p,model.indexlevel{p},:) = 1; %(eq.1 ) to incorporate local level
                
                
                if settings.cycle
                    for i=1:settings.numbercycles    
                        Z(p,model.indexcycles{p}(i*2-1),:) = 1; % (eq.1 and eq.4) to incorporate all cycle terms c_t (c1_t...cn_t)
                        %   Z(:,model.indexcycles{p}(1)+(i-1)*2) = 1; % (eq.1 and eq.4) to incorporate all cycle terms c_t (c1_t...cn_t)
                    end
                end
                
                if settings.regression
                    for i=1:settings.numberregressors
                        
                        Z(p,model.indexregressors{p}(i),:)=tseries.regressors(i,:);
                        
                    end
                end
                
                if settings.intervention
                    for i=1:settings.numberinterventions
                        Z(p,model.indexinterventions{p}(i),1:settings.tau(i)-1)=0;
                        Z(p,model.indexinterventions{p}(i),settings.tau(i):end)=1;
                    end
                end
                
                if settings.AR
                    Z(p,model.indexAR{p},:) = 1;
                    if settings.ARorder>1
                        for i=2:settings.ARorder
                            Z(p,model.indexAR{p}+i-1,:) = 0;
                        end
                    end
                end
                
                % univariate
                % including quadratic and two cycle terms Z should be: [1 0 1 0 0 1 0 1 0]
                
                % multivariate
                % including quadratic and two cycle terms Z should be: [1 0 1 0 0 1 0 1 0    0 0 0 0 0 0 0 0 0;
                %                                                       0 0 0 0 0 0 0 0 0    1 0 1 0 0 1 0 1 0]
                
            end
            
            %% T: transition matrix
            
            % loop over time
            for i=1:tseries.ntimes-1
                % loop over number of time series
                for p=1:tseries.ntseries
                    
                    if settings.acc
                        % quadratic term (eq.3)
                        T(model.indexlevel{p}:model.indexacc{p},model.indexlevel{p}:model.indexacc{p},i)=  [1  dt(i)  0.5*dt(i)^2;
                            0  1  dt(i);
                            0  0  1];
                    else
                        if settings.slope
                            % trend (eq.2)
                            T(model.indexlevel{p}:model.indexslope{p},model.indexlevel{p}:model.indexslope{p},i)=  [1 dt(i);
                                0 1];
                        else
                            % local level
                            T(model.indexlevel{p},model.indexlevel{p})=1;
                        end
                    end
                    
                    % cycle terms (eq.4)
                    if settings.cycle
                        for ii=1:settings.numbercycles
                            T(model.indexcycles{p}(ii*2-1):model.indexcycles{p}(ii*2),...
                                model.indexcycles{p}(ii*2-1):model.indexcycles{p}(ii*2),i)= ...
                                [cos(settings.lambda(ii)*dt(i)) sin(settings.lambda(ii)*dt(i));
                                -sin(settings.lambda(ii)*dt(i)) cos(settings.lambda(ii)*dt(i))];
                        end
                    end
                    
                    if settings.regression
                        for ii=1:settings.numberregressors
                            T(model.indexregressors{p}(ii),model.indexregressors{p}(ii),i)=1;
                        end
                    end
                    
                    if settings.intervention
                        for ii=1:settings.numberinterventions
                            T(model.indexinterventions{p}(ii),model.indexinterventions{p}(ii),i)=1;
                        end
                    end
                    
                    if settings.AR
                        if isempty(settings.ARorder) || settings.ARorder < 0
                            
                            error(strcat('settings.ARorder has an invalid value:',num2str(settings.ARorder)))
                        end
                        if isempty(settings.ARphi)
                            error('set settings.ARphi to a valid value')
                        end
                        if settings.ARorder>1
                            t_input=zeros(settings.ARorder,1);
                            t_input(2) = 1;
                            T_ar = toeplitz(zeros(settings.ARorder,1),t_input);
                            for ii=1:settings.ARorder
                                T_ar(ii,1) = settings.ARphi(p,ii);
                            end
                            T(model.indexAR{p}:model.indexAR{p}+settings.ARorder-1,model.indexAR:model.indexAR+settings.ARorder-1) = T_ar;
                        else
                            if numel((settings.ARphi)) == 1
                                T(model.indexAR{p},model.indexAR{p},i)=settings.ARphi;
                            else
                                T(model.indexAR{p},model.indexAR{p},i)=settings.ARphi(p,1);
                            end
                        end
                    end
                end
            end
            
            %% Q: process noise
            
            
            
            %settings.processnoise='Diagonal';
            for i=1:tseries.ntimes-1
                % loop over number of time series
                for p=1:tseries.ntseries
                    
                    % level
                    Q(model.indexlevel{p},model.indexlevel{p},i)=distvar.level(p)*dt(i);
                    % define full process noise variance vector (without
                    % time dependence)
                    if (i==1)
                        distvar.eta(model.indexlevel{p})=distvar.level(p);
                    end
                    
                    if settings.multivariate
                        for pp=p+1:tseries.ntseries
                            Q(model.indexlevel{p},model.indexlevel{pp},i)=distvar.covlevel(p,pp)*dt(i);
                            Q(model.indexlevel{pp},model.indexlevel{p},i)=distvar.covlevel(p,pp)*dt(i);
                        end
                    end
                    
                    if settings.slope
                        % trend
                        Q(model.indexslope{p},model.indexslope{p},i)= distvar.slope(p)*dt(i);
                        % define full process noise variance vector (without
                        % time dependence)
                        if (i==1)
                            distvar.eta(model.indexslope{p})=distvar.slope(p);
                        end
                        
                        if settings.multivariate
                            for pp=p+1:tseries.ntseries
                                Q(model.indexslope{pp},model.indexslope{p},i)=distvar.covslope(p,pp)*dt(i);
                                Q(model.indexslope{p},model.indexslope{pp},i)=distvar.covslope(p,pp)*dt(i);
                            end
                        end
                    end
                    
                    if settings.acc
                        % quadratic term
                        Q(model.indexacc{p},model.indexacc{p},i)=distvar.acc(p)*dt(i);
                        % define full process noise variance vector (without
                        % time dependence)
                        if (i==1)
                            distvar.eta(model.indexslope{p})=distvar.slope(p);
                        end
                        
                        if settings.multivariate
                            for pp=p+1:tseries.ntseries
                                Q(model.indexacc{p},model.indexacc{pp},i)=distvar.covacc(p,pp)*dt(i);
                                Q(model.indexacc{pp},model.indexacc{p},i)=distvar.covacc(p,pp)*dt(i);
                            end
                        end
                    end
                    
%                     % process noise regressor
%                     if settings.regression
%                         if settings.timevariableregressor
%                             for ii=1:settings.numberregressors
%                                 Q(model.indexregressors{p},model.indexregressors{p},i)=distvar.regressor(p,ii)*dt(i);
%                                 if (i==1)
%                                     distvar.eta(model.indexregressors{p})=distvar.regressor(p,ii);
%                                 end
%                             end
%                         end
%                     end
%                     
                    
                    
                    
                    % cycle terms, time dependency see Harvey eq. 9.2.7
                    if settings.cycle
                        
                        % use different variances for omega_t and omega_t*
                        for ii=1:settings.numbercycles
                            Q(model.indexcycles{p}(ii*2-1),model.indexcycles{p}(ii*2-1),i) = distvar.cycle(p,ii)*dt(i);
                            Q(model.indexcycles{p}(ii*2),model.indexcycles{p}(ii*2),i) = distvar.cyclestar(p,ii)*dt(i);
                            Q(model.indexcycles{p}(ii*2-1),model.indexcycles{p}(ii*2),i)= distvar.covcyclestar(p,ii)*dt(i);
                            Q(model.indexcycles{p}(ii*2),model.indexcycles{p}(ii*2-1),i)= distvar.covcyclestar(p,ii)*dt(i);
                            
                            % define full process noise variance vector (without
                            % time dependence)
                            if (i==1)
                                distvar.eta(model.indexcycles{p}(ii*2-1))=distvar.cycle(p,ii);
                                distvar.eta(model.indexcycles{p}(ii*2))=distvar.cyclestar(p,ii);
                            end
                        end

                    end
                    
                    % autoregressive term
                    if settings.AR
%                         if isempty(distvar.AR)
%                             error('provide initial process noise for AR')
%                         end
                        for ii=1:settings.ARorder
                           
                            Q(model.indexAR{p}+ii-1,model.indexAR{p}+ii-1,i) = distvar.AR(p)*dt(i);
                            distvar.eta(model.indexAR{p}+ii-1) = distvar.AR(p);
                        end
                        
                        if settings.multivariate
                            for pp=p+1:tseries.ntseries
                                Q(model.indexAR{p},model.indexAR{pp},i)=distvar.covAR(p,pp)*dt(i);
                                Q(model.indexAR{pp},model.indexAR{p},i)=distvar.covAR(p,pp)*dt(i);
                            end
                        end
                    end
                end
            end
            
            % standard definition
            R=eye(n);
            
            if settings.AR
                if settings.ARorder>1
                    for ii=1:settings.ARorder-1
                        R(model.indexAR+ii,model.indexAR+ii) = 0;
                    end
                end
            end
            
            % full variance covariance matrix for process noise 
            % without time dependence
            distvar.varcoveta=Q(:,:,1)/tseries.dt(1);
            
        else
            error('invalid choice for state space formulation')
        end
        
        if settings.multivariate
           % make matrix Wt to filter out missing observations
           % see 4.10 of Durbin and Koopman
           
             model.W=~isnan(tseries.y)'; 
         % unsure    % we multiply W(i,:) with I to create Wt as described in
             % Durbin and Koopman when Wt is needed
        end

        
    case 0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% integer time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % case with integer time steps
        
        if settings.regression
            error('regression not implemented for integer time')
        end
        if settings.intervention
            error('intervenion not implemented for integer time')
        end
        if settings.AR
            error('autoregression not implemented for integer time')
        end
        if settings.multivariate
            error('multivariateness not implemented for integer time')
        end
        
        if isempty(settings.type);settings.type='Normal_Slope_and_Acc';end
        
        
        
        %     if strcmp(settings.type,'Lag_Operator_Both_Slope_and_Acc')
        %         %% Determine dimension of matrices
        %         % trend adds 2 dimensions
        %         % acceleration 3
        %         % each cycle 2
        %         if ~settings.slope
        %             disp('When lag operator is used settings.slope has to be true')
        %             return
        %         end
        %
        %         n=2;
        %
        %         if settings.acc
        %             n=n+3;
        %         end
        %         if settings.cycle
        %             n=n+2*settings.numbercycles;
        %         end
        %
        %         %% make matrices
        %         Z = zeros(1,n);
        %         T = zeros(n,n);
        %         Q = zeros(n,n);
        %
        %         %% H: disturbance variance
        %         H = eye(settings.ntimes,settings.ntimes)*distvar.irr;
        %
        %         %% Z: design matrix
        %
        %         Z(1) = 1; %(eq.1 and eq.2) to incorporate trend mu_t+1
        %
        %         NextDim=3;
        %
        %         if settings.acc
        %             Z(NextDim) = 1; %(eq.1 and eq.3) to incorporate quadratic term mu2_t+1
        %             NextDim=NextDim+3;
        %         end
        %
        %
        %         if settings.cycle
        %            CycleDim=NextDim;
        %            for i=1:settings.numbercycles
        %               NextDim=NextDim+2;
        %         %      nStartCycle(i)=n+1-(settings.numbercycles-i+1)*2;
        %               Z(CycleDim+(i-1)*2) = 1; % (eq.1 and eq.4) to incorporate all cycle terms c_t (c1_t...cn_t)
        %            end
        %         end
        %
        %         % including quadratic and two cycle terms Z should be: [1 0 1 0 0 1 0 1 0]
        %
        %         %% T: transition matrix
        %
        %         % trend (eq.2)
        %         T(1:2,1:2)=[2 -1;
        %                     1 0];
        %         % quadratic term (eq.3)
        %         if settings.acc
        %            T(3:5,3:5)= [3 -3 1;
        %                         1  0 0;
        %                         0  1 0];
        %         end
        %
        %         % cycle terms (eq.4)
        %         if settings.cycle
        %            for i=1:settings.numbercycles
        %                T(CycleDim+(i-1)*2:CycleDim+(i-1)*2+1,CycleDim+(i-1)*2:CycleDim+(i-1)*2+1)=[cos(settings.lambda(i)) sin(settings.lambda(i));
        %                                                                                   -sin(settings.lambda(i))  cos(settings.lambda(i))];
        %            end
        %         end
        %
        %
        %         %% Q: process noise
        %
        %         % trend
        %         Q(1,1)=distvar.slope;
        %         % define full process noise variance vector
        %         distvar.eta(1)=distvar.slope;
        %         % quadratic term
        %         if settings.acc
        %            Q(3,3)=distvar.acc;
        %            distvar.eta(3)=distvar.acc
        %         end
        %
        %         % cycle terms
        %         if settings.cycle
        %             for i=1:settings.numbercycles
        %                 Q(CycleDim+(i-1)*2,CycleDim+(i-1)*2)=distvar.cycle(i);
        %                 Q(CycleDim+(i-1)*2+1,CycleDim+(i-1)*2+1)=distvar.cycle(i);
        %             end
        %         end
        %
        %
        %         R=eye(n);
        %
        %     elseif strcmp(settings.type,'Lag_Operator_either_Slope_or_Acc')
        %         %% Determine dimension of matrices
        %         % trend uses 2 dimensions
        %         % or instead acceleration 3
        %         % each cycle 2
        %         if ~settings.slope
        %             disp('When lag operator is used settings.slope has to be true')
        %             return
        %         end
        %
        %
        %
        %         if settings.acc
        %             n=3;
        %         else
        %             n=2;
        %         end
        %         if settings.cycle
        %             n=n+2*settings.numbercycles;
        %         end
        %
        %         %% make matrices
        %         Z = zeros(1,n);
        %         T = zeros(n,n);
        %         Q = zeros(n,n);
        %
        %         %% H: disturbance variance
        %         H = eye(settings.ntimes,ntimes)*distvar.irr;
        %
        %         %% Z: design matrix
        %
        %
        %
        %         if settings.acc
        %             % acceleration only has random walk
        %             Z(1) = 1; %(eq.1 and eq.3) to incorporate quadratic term mu2_t+1
        %             NextDim=4;
        %         else
        %             % slope only has random walk
        %             Z(1) = 1; %(eq.1 and eq.2) to incorporate trend mu_t+1
        %             NextDim=3;
        %         end
        %
        %         if settings.cycle
        %            CycleDim=NextDim;
        %            for i=1:settings.numbercycles
        %               NextDim=NextDim+2;
        %         %      nStartCycle(i)=n+1-(settings.numbercycles-i+1)*2;
        %               Z(CycleDim+(i-1)*2) = 1; % (eq.1 and eq.4) to incorporate all cycle terms c_t (c1_t...cn_t)
        %            end
        %         end
        %
        %         % including quadratic and two cycle terms Z should be: [1 0 1 0 0 1 0 1 0]
        %
        %
        %         %% T: transition matrix
        %
        %
        %         if settings.acc
        %         % quadratic term (eq.3)
        %            T(1:3,1:3)= [3 -3 1;
        %                         1  0 0;
        %                         0  1 0];
        %         else
        %         % trend (eq.2)
        %             T(1:2,1:2)=[2 -1;
        %                         1 0];
        %         end
        %
        %         % cycle terms (eq.4)
        %         if settings.cycle
        %            for i=1:settings.numbercycles
        %                T(CycleDim+(i-1)*2:CycleDim+(i-1)*2+1,CycleDim+(i-1)*2:CycleDim+(i-1)*2+1)=[cos(settings.lambda(i)) sin(settings.lambda(i));
        %                                                                                   -sin(settings.lambda(i))  cos(settings.lambda(i))];
        %            end
        %         end
        %
        %         %% Q: process noise
        %
        %
        %
        %         if settings.acc
        %             % quadratic term
        %            Q(1,1)=distvar.acc;
        %         else
        %             % trend
        %             Q(1,1)=distvar.slope;
        %         end
        %
        %         % cycle terms
        %         if settings.cycle
        %             for i=1:settings.numbercycles
        %                 Q(CycleDim+(i-1)*2,CycleDim+(i-1)*2)=distvar.cycle(i);
        %                 Q(CycleDim+(i-1)*2+1,CycleDim+(i-1)*2+1)=distvar.cycle(i);
        %             end
        %         end
        %
        % %         disp('Process Noise')
        % %         Q
        %
        %         R=eye(n);
        
        if strcmp(settings.type,'Normal_Slope_and_Acc')
            %% Determine dimension of matrices
            % trend uses 2 dimensions (mu + v)
            % acceleration adds 1 extra
            % each cycle 2
            
            %% Determine dimension of matrices
            model.indexslope=[];
            model.indexacc=[];
            model.indexcycles=[];
            model.indexregressors=[];
            model.indexinterventions=[];
            model.indexAR=[];
            % trend uses 2 dimensions (mu + v)
            % acceleration adds 1 extra
            % each cycle 2
            % each regressor 1
            % each intervention 1
            
            % loop on number of time series
            n=0;
            
            n=n+1; % current dimension count
            model.indexlevel{1} = n;
            
            if settings.slope
                n=n+1;
                model.indexslope{1} = n;
                
                if settings.acc
                    n=n+1;
                    model.indexacc{1} = n;
                end
            end
            
            if settings.cycle
                model.indexcycles{1}=[n+1:n+2*settings.numbercycles];
                %   model.indexcycles=[model.indexcycles n+1:n+2*settings.numbercycles];
                n=n+2*settings.numbercycles;
            end
            
          
          
            
            %% make matrices
            Z = zeros(1,n);
            T = zeros(n,n);
            Q = zeros(n,n);
            
            %% H: disturbance variance
            H = eye(tseries.ntimes,tseries.ntimes)*distvar.irr;
            
            %% Z: design matrix
            Z(1) = 1; %(eq.1 ) to incorporate local level
            %         if settings.slope
            %             NextDim=3;
            %             if settings.acc
            %                 NextDim=4;
            %             end
            %         else
            %             % local level model
            %             NextDim=2;
            %         end
            
            if settings.cycle
                
                for i=1:settings.numbercycles
                    %  NextDim=NextDim+2;
                    Z(1,model.indexcycles{1}+(i-1)*2) = 1; % (eq.1 and eq.4) to incorporate all cycle terms c_t (c1_t...cn_t)
                end
            end
            
            %         if settings.cycle
            %            CycleDim=NextDim;
            %            for i=1:settings.numbercycles
            %               NextDim=NextDim+2;
            %         %      nStartCycle(i)=n+1-(settings.numbercycles-i+1)*2;
            %               Z(CycleDim+(i-1)*2) = 1; % (eq.1 and eq.4) to incorporate all cycle terms c_t (c1_t...cn_t)
            %            end
            %         end
            
            % including quadratic and two cycle terms Z should be: [1 0 1 0 0 1 0 1 0]
            
            %% T: transition matrix
            
            
            if settings.acc
                % quadratic term (eq.3)
                T(1:3,1:3)= [1  1  0.5;
                    0  1  1;
                    0  0  1];
            else
                if settings.slope
                    % trend (eq.2)
                    T(1:2,1:2)=[1 1;
                        0 1];
                else
                    % local level
                    T(1,1)=1;
                end
            end
            
            % cycle terms (eq.4)
            if settings.cycle
                for i=1:settings.numbercycles
                    T(model.indexcycles{1}+(i-1)*2:model.indexcycles{1}+(i-1)*2+1,model.indexcycles{1}+(i-1)*2:model.indexcycles{1}+(i-1)*2+1)=[cos(settings.lambda(i)) sin(settings.lambda(i));
                        -sin(settings.lambda(i))  cos(settings.lambda(i))];
                end
            end
            
            %% Q: process noise
            
            settings.processnoise='Diagonal';
            
            if settings.acc
                if strcmp(settings.processnoise,'Full')
                    disp('full process noise not implemented when accelerations are included')
                else
                    % quadratic, slope and level variance terms
                    Q(3,3)=distvar.acc;
                    Q(2,2)=distvar.slope;
                    Q(1,1)=distvar.level;
                    distvar.eta(3)=distvar.acc;
                    distvar.eta(2)=distvar.slope;
                    distvar.eta(1)=distvar.level;
                end
            else
                Q(1,1)=distvar.level;
                distvar.eta(1)=distvar.level;
                if settings.slope
                    % trend and level terms
                    Q(2,2)=distvar.slope;
                    distvar.eta(2)=distvar.slope;
                    if strcmp(settings.processnoise,'Full')
                        Q(1,1)=distvar.slope*1/3+distvar.level;
                        Q(2,1)=distvar.slope*0.5;
                        Q(1,2)=distvar.slope*0.5;
                    end
                end
            end
            
            if settings.cycle
                % check the dimension of distvar.cycle
                [CycleRow CycleCol]=size(distvar.cycle);
                % cycle terms
                if CycleRow==1
                    % use identical variances for omega_t and omega_t*
                    for i=1:settings.numbercycles
                           Q(model.indexcycles{1}(i*2-1),model.indexcycles{1}(i*2-1),i) = distvar.cycle(1,i);
                            Q(model.indexcycles{1}(i*2),model.indexcycles{1}(i*2),i) = distvar.cycle(1,i);
                         %   Q(model.indexcycles{p}(ii*2-1),model.indexcycles{p}(ii*2),i)= distvar.covcyclestar(p,ii)*dt(i);
                          %  Q(model.indexcycles{p}(ii*2),model.indexcycles{p}(ii*2-1),i)= distvar.covcyclestar(p,ii)*dt(i);
%                         Q(model.indexcycles{1}+(i-1)*2,model.indexcycles{1}+(i-1)*2)=distvar.cycle(i);
%                         Q(model.indexcycles{1}+(i-1)*2+1,model.indexcycles{1}+(i-1)*2+1)=distvar.cycle(i);
                    end
                elseif CycleRow==2
                    % use different variances for omega_t and omega_t*
                    for i=1:settings.numbercycles
                        Q(model.indexcycles{1}(i*2-1),model.indexcycles{1}(i*2-1),i) = distvar.cycle(1,i);
                            Q(model.indexcycles{1}(i*2),model.indexcycles{1}(i*2),i) = distvar.cyclestar(1,i);
%                         Q(model.indexcycles{1}+(i-1)*2,model.indexcycles{1}+(i-1)*2)=distvar.cycle(1,i);
%                         Q(model.indexcycles{1}+(i-1)*2+1,model.indexcycles{1}+(i-1)*2+1)=distvar.cycle(2,i);
                    end
%                     

%                     
                else
                    disp('distvar.cycle has too many rows, exit')
                    return
                    
                end
            end
            % standard definition
            R=eye(n);
            
            % full variance covariance matrix for process noise 
            % without time dependence
            distvar.varcoveta=Q;
        else
            disp('invalid choice for state space formulation')
            settings.type
            return
        end
        
end

% save in structure array
model.dim = n;
model.Z = Z;
model.T = T;
model.Q = Q;
model.H = H;
model.R = R;

end