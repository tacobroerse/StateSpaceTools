function [estfilter] = Kalman(tseries,model,settings,estfilter)
% Kalman.m performs Kalman filtering. The Kalman filter serves as a update of our knowledge of the
% state each time a new observation y_t is brought in. Time dependent
% matrix version.

%
% HOW:     [estfilter] = Kalman(tseries,model,settings,estfilter)
%
% Input:
%       tseries, structure of which this function uses:
%           tseries.t [epochs x 1]          normalised time vector
%           tseries.Y                       observation vector
%           tseries.ntimes                  number of epochs
%           tseries.missingepochs           indexes of missing data (only for integer time)
%
%       model, structure of which this function uses:
%           model.Z   [n x epochs]      design matrix
%           model.T   [n x n x epochs]      transition matrix
%           model.Q   [n x n x epochs]      process noise matrix
%           model.R   [n x n]               selection matrix
%           model.H   [epochs]              disturbance variance matrix
%           model.dim                       dimension of the state vector
%
%       settings, structure of which this function uses
%           settings.filter                 type of Kalman formulation
%           settings.continuoustime         whether time is non-integer (thus continuous)
%
%       estfilter, structure of which this function uses
%           estfilter.a [n x t]             the one step ahead prediction of alpha, initialized for t=1: a(:,1)
%           estfilter.P [n x n x t]         variance of a, initialized for t=1: P(:,:,1)
%
% Output:
%       estfilter, structure of which this function outputs
%           estfilter.a [n x epochs]        filtered state
%           estfilter.P [n x n  x epochs]   state error variance
%           estfilter.K [n x epochs]        Kalman gain
%           estfilter.v [epochs]            Prediction error
%           estfilter.F [epochs]            The variance of prediction error v
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
% uses: householder.m
%
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 June 2014 DBT Broerse
%
% Version 1.1 September 2014 DBT Broerse
% - design matrix now time dependent
%
% Version 1.2 October 2014 Thomas Frederikse
% - added: square root filter
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
% - continuous time and integer time functions merged
%
% Version 2.1 March 2015 DBT Broerse
% - changed size of Z matrix
%
% Version 2.2 June 2018 DBT Broerse
% - changed order of Z matrix for continuous time
% - adapted the filter for multivariate models
%
% Version 2.3 May 2019 DBT Broerse
% - addition of missing epochs in multivariate models
%
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------
%
%%%% THEORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Formulation of the Kalman filter according to chapter 4.3.2 of Durbin and
% Koopman, equation 4.24.
%
% Kalman filtering consist of following recursion for t=1,2..,N
%
% Difference between observation and one step ahead prediction
% v_t = y_t - Z_t a_t                                               (eq.1)
%
% where v_t is the prediction error, y_t is the observation at time t,
% Z is the design matrix, a_t is the one step ahead prediction of alpha,
% the state vector.
%
% The variance of prediction error v_t
% F_t = Z_t P_t Z_t' + H_t                                          (eq.2)
%
% where H_t is the observational noise
%
% Filtered estimator can be skipped:
% Determination of filtered estimator at_t (at|t in Durbin & Koopman)
% at_t = a_t + P_t Z_t' F_t^-1 v_t                                 (eq.3)
%
% Determination of filtered estimator variance Pt_t (Pt|t in Durbin & Koopman)
% Pt_t = P_t - P_t Z_t' F_t^-1 Z_t P_t                             (eq.4)
%
% Proceed here again:
% Kalman gain K
% K_t = T_t P_t Z_t' inv(F_t)                                       (eq.5)
%
% Prediction step
% a_t+1 = T_t a_t + K_t v_t                                         (eq.6)
%
% P_t+1 = T_t P_t (T_t - K_t Z_t)' + R_t Q_t R_t'                           (eq.7)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default
if ~any(strcmp('filter',fieldnames(settings)));settings.filter=[];end
if isempty(settings.filter);settings.filter='DurbinKoopman';end

% extract from structure
a=estfilter.a;
P=estfilter.P;
y=tseries.Y;

Z=model.Z;
Q=model.Q;
R=model.R;
T=model.T;
H=model.H;

if settings.multivariate
    % mask matrix for multivariate timeseries with missing observations
    W=model.W;
end

% different code for continuous time (general case) and integer time steps
switch settings.continuoustime
    
    case 1
        % continuous time
        %
        % initialise
        if settings.multivariate
            v=zeros(tseries.ntseries,tseries.ntimes);
            F=zeros(tseries.ntseries,tseries.ntseries,tseries.ntimes);
            K=zeros(model.dim,tseries.ntseries,tseries.ntimes);
        else
            v=zeros(tseries.ntimes,1);
            F=zeros(tseries.ntimes,1);
            K=zeros(model.dim,tseries.ntimes);
        end
        
        if strcmp(settings.filter,'DurbinKoopman')
            for i=1:tseries.ntimes
                
                if settings.multivariate
                    % multivariate case
                    Zt=Z(:,:,i);
                    Ht=H(:,:,i);
                    yt=y(:,i);
                    
                    
                    if tseries.epochsmissingobs(i)
                        % if there are any missing observations
                        % special treatment for epochs with missing observations:
                        % set Z to 0 (page 111 of Durbin & Koopman)
                        % we take a slightly different approach than what
                        % is written at page 111, instead of masking out
                        % the dimensions of y, H and Z that are not
                        % observed, we remove those dimensions, which is
                        % what Durbin and Koopman actually seem to suggest
             
                        Windex=find(W(i,:));
                     
                        Ztstar=Zt(Windex,:);
                        Htstar=Ht(Windex,Windex);
                        ytstar=yt(Windex);
%                         Ztstar=Wt*Zt;% mask out dimensions that are not observed
%                         Htstar=Wt*Ht*Wt';% mask out dimensions that are not observed
%                         ytstar=yt;
%                         ytstar(~Windex)=0; % replace NaN
                    
                        
                        v(Windex,i)=ytstar-Ztstar*a(:,i); % (eq.1)
                        % with variance
                        F(Windex,Windex,i)=Ztstar*P(:,:,i)*Ztstar'+Htstar; % (eq.2)
                        
                        % Kalman gain (if Wt contains zeros, only update Kalman
                        % gain for time series with observations)
                        % reduce dimension of Zt and F here, to keep only
                        % indexes with observations, otherwise
                        % singularities will arise from inverse of F
                     
                        % only compute K for observed dimensions
                        K(:,Windex,i)=T(:,:,i)*P(:,:,i)*Ztstar'/(F(Windex,Windex,i)); % (eq.5)
                        
                    else
                        % difference between observation and model
                        %(a.k.a. one step ahead prediction error)
                        
                        v(:,i)=yt-Zt*a(:,i); % (eq.1)
                        % with variance
                        F(:,:,i)=Zt*P(:,:,i)*Zt'+Ht; % (eq.2)
                        % Kalman gain 
                        K(:,:,i)=T(:,:,i)*P(:,:,i)*Zt'/(F(:,:,i)); % (eq.5)
                    end
                    
                    
                    
                    % prediction
                    if  i<tseries.ntimes
                        a(:,i+1)=T(:,:,i)*a(:,i)+K(:,:,i)*v(:,i);
                        P(:,:,i+1)=T(:,:,i)*P(:,:,i)*(T(:,:,i)-K(:,:,i)*Zt)'+R * Q(:,:,i)*R';
                    end
                    
                else
                    
                    % univariate case
                    Zt=Z(1,:,i);
                    Ht=H(i,i);
                    yt=y(:,i);
                    
                    if ~isempty(find(i==tseries.missingepochs))
                        % special treatment for epochs with missing observations:
                        % set Z to 0 (page 111 of Durbin & Koopman)
                        % this reduces v, K and F to 0
                        % prediction
                        if  i<tseries.ntimes
                            a(:,i+1)=T(:,:,i)*a(:,i);
                            P(:,:,i+1)=T(:,:,i)*P(:,:,i)*T(:,:,i)'+R*Q(:,:,i)*R';
                        end
                    else
                        % difference between observation and model
                        % (a.k.a. one step ahead prediction error)
                        v(i)=yt-Zt*a(:,i); % (eq.1)
                        % with variance
                        F(i)=Zt*P(:,:,i)*Zt'+Ht; % (eq.2)
                        % Kalman gain
                        K(:,i)=T(:,:,i)*P(:,:,i)*Zt'/(F(i)); % (eq.5)
                        
                        
                        
                        % prediction
                        if  i<tseries.ntimes
                            a(:,i+1)=T(:,:,i)*a(:,i)+K(:,i)*v(i);
                            P(:,:,i+1)=T(:,:,i)*P(:,:,i)*(T(:,:,i)-K(:,i)*Zt)'+R * Q(:,:,i)*R';
                        end
                    end
                end
                
            end
        elseif strcmp(settings.filter,'SquareRoot')
            % square root filter to avoid negative variances due to numeric
            % errors, see Durbin & Koopman chapter 6.3
            for i=1:tseries.ntimes
                if settings.multivariate
                    if i==1
                        disp('untested combination: square root filter and multivariate models')
                    end
                    Zt=Z(:,:,i);
                    Ht=H(:,:,i);
                    yt=y(:,i);
                    
                    %(one step ahead prediction error)
                    v(:,i)=yt-Zt*a(:,i);
                    % produce a lower triangular matrix
                    Ptilde(:,:) = chol(P(:,:,i),'lower');
                    Qtilde(:,:) = sqrt(Q(:,:,i));
                    Htilde      = sqrt(Ht);
                    Uparta = [Zt*Ptilde(:,:)  Htilde             zeros(tseries.ntseries, model.dim)];
                    Upartb = [T(:,:,i)*Ptilde(:,:)  zeros(model.dim,tseries.ntseries)  R(:,:)*Qtilde(:,:)];
                    U      = [Uparta;Upartb];
                    Ustar  = householder(U(:,:));
                    ustar1 = Ustar(1,1);
                    ustar2 = Ustar(2:model.dim+1,1);
                    ustar3 = Ustar(2:model.dim+1,2:model.dim+1);
                    Ptp1   = ustar3*ustar3';
                    atp1 = T(:,:,i)*a(:,i) + ustar2/ustar1*v(i);
                    Ft   = ustar1*ustar1';
                    Kt   = T(:,:,i)*P(:,:,i)*Zt'/(Ft);
                    F(:,:,i) = Ft;
                    % Kalman gain
                    K(:,:,i) = Kt;
                    if  i<tseries.ntimes
                        a(:,i+1)   = atp1;
                        P(:,:,i+1) = Ptp1;
                    end
                else
                    Zt=Z(1,:,i);
                    Ht=H(i,i);
                    yt=y(i);
                    
                    %(one step ahead prediction error)
                    v(i)=yt-Zt*a(:,i);
                    % produce a lower triangular matrix
                    Ptilde(:,:) = chol(P(:,:,i),'lower');
                    Qtilde(:,:) = sqrt(Q(:,:,i));
                    Htilde      = sqrt(Ht);
                    Uparta = [Zt*Ptilde(:,:)  Htilde             zeros(1, model.dim)];
                    Upartb = [T(:,:,i)*Ptilde(:,:)  zeros(model.dim,1)  R(:,:)*Qtilde(:,:)];
                    U      = [Uparta;Upartb];
                    Ustar  = householder(U(:,:));
                    ustar1 = Ustar(1,1);
                    ustar2 = Ustar(2:model.dim+1,1);
                    ustar3 = Ustar(2:model.dim+1,2:model.dim+1);
                    Ptp1   = ustar3*ustar3';
                    atp1 = T(:,:,i)*a(:,i) + ustar2/ustar1*v(i);
                    Ft   = ustar1*ustar1';
                    Kt   = T(:,:,i)*P(:,:,i)*Zt'/(Ft);
                    F(i) = Ft;
                    % Kalman gain
                    K(:,i) = Kt;
                    if  i<tseries.ntimes
                        a(:,i+1)   = atp1;
                        P(:,:,i+1) = Ptp1;
                    end
                end
            end
        else
            disp('invalid filter type')
            return
        end
        
    case 0
        
        % integer time steps
        
        % initialise
        v=zeros(tseries.ntimes,1);
        F=zeros(tseries.ntimes,1);
        K=zeros(model.dim,tseries.ntimes);
        
        if strcmp(settings.filter,'DurbinKoopman')
            for i=1:tseries.ntimes
                % check if i is a missing epoch
                if ~isempty(find(i==tseries.missingepochs))
                    % special treatment for epochs with missing observations:
                    % set Z to 0 (page 111 of Durbin & Koopman)
                    if  i<tseries.ntimes
                        a(:,i+1)=T*a(:,i);
                        P(:,:,i+1)=T*P(:,:,i)*T'+R*Q*R';
                    end
                    
                else
                    % difference between observation and model
                    %(a.k.a. one step ahead prediction error)
                    
                    v(i)=y(i)-Z*a(:,i); % (eq.1)
                    % with variance
                    F(i)=Z*P(:,:,i)*Z'+H(i,i); % (eq.2)
                    % Kalman gain
                    K(:,i)=T*P(:,:,i)*Z'/(F(i)); % (eq.5)
                    %K(:,i)=T*P(:,:,i)*Z'*inv(F(i)); % (eq.5)
                    % prediction
                    if  i<tseries.ntimes
                    
                        a(:,i+1)=T*a(:,i)+K(:,i)*v(i);
                        P(:,:,i+1)=T*P(:,:,i)*(T-K(:,i)*Z)'+R*Q*R';
                    end
                    
                end
                
            end
        else
            error('only filter DurbinKoopman is implemented for integer time')
        end
end

% save in structure
estfilter.a=a;
estfilter.P=P;
estfilter.K=K;
estfilter.v=v;
estfilter.F=F;

end

