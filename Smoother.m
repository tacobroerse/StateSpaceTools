function [estsmooth] = Smoother(tseries,settings,model,estfilter)
% SmootherContTime.m performs Kalman filtering. The smoother recursively 
% smoothes the Kalman derived estimates using all available observations. 
% General version, for continuous time (any time step).

%
% HOW:     [estsmooth] = Smoother(tseries,settings,model,estfilter)
%
% Input:  
%       tseries, structure of which this function uses:
%           tseries.ntimes                  number of epochs  
%
%       model, structure of which this function uses:
%           model.Z   [epochs x n]      design matrix
%           model.T   [n x n x epochs]      transition matrix
%           model.dim                       dimension of the state vector
%
%       settings, structure of which this function uses
%           settings.smoother               type of smoother formulation
%           settings.continuoustime         whether time is non-integer (thus continuous)
%           settings.multivariate
%
%       estfilter, structure of which this function uses
%           estfilter.a   [n x t]           the one step ahead prediction of alpha
%           estfilter.P   [n x n x t]       variance of a
%           estfilter.K [n x epochs]        Kalman gain
%           estfilter.v [epochs]            Prediction error
%           estfilter.F [epochs]            The variance of prediction error v
%
%   
% Output:   
%       estsmooth, structure of smoothed estimates    
%           estsmooth.alpha [n x t]         smoothed state vector
%           estsmooth.V [n x n x t]         smoothed state error variance                    
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
% - added: square root smoother
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
% - continuous time and integer time functions merged
%
% Version 2.1 June 2018 DBT Broerse
% - extended to multivariate models
% - changed indexing of Z for continuous time
%
% Version 2.2 June 2019 DBT Broerse
% - addition of missing epochs in multivariate models
%
%----------------------------------------------------------------------------
% remarks: 
%----------------------------------------------------------------------------
%
%%%% THEORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Formulation of the state smoothing according to chapter 4.4.4 of Durbin and
% Koopman, equation 4.44. 
%
% Smoothing consist of following recursion for t=n,..,1
%
% Difference between observation and one step ahead prediction 
% v_t = y_t - Z_t a_t                                               (eq.1)
%
% The weighted sum of innovations r, with r_n = 0
% r_t-1 = Z_t' inv(F_t) v_t + L_t' r_t                              (eq.2)
%
% Weighted sum of the inverse variances of innovations N, with N_n = 0
%
% N_t-1 = Z_t' inv(F_t) Z_t + L_t' N_t L_t                          (eq.3)
%
% where L_t = T_t - K_t Z_t                                         (eq.4)
% and K_t is the Kalman gain, Z_t the design
% matrix and T_t the transition matrix
%
% The smoothed state vector alpha becomes:
% alpha_t = a_t + P_t r_t-1                                         (eq.5)
%
% and the smoothed state variance V_t:
% V_t = P_t - P_t N_t-1 P_t
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default
if ~any(strcmp('smoother',fieldnames(settings)))
    settings.smoother=[];
end
if isempty(settings.smoother)
    if isempty(settings.filter)
        settings.smoother='DurbinKoopman';
    else
        settings.smoother=settings.filter;
    end
end


% extract from structure
a=estfilter.a;
P=estfilter.P;
K=estfilter.K;
v=estfilter.v;
F=estfilter.F;

Z=model.Z;
T=model.T;
if settings.multivariate
    % mask for missing observations
    W=model.W;
end


% different code for continuous time (general case) and integer time steps
switch settings.continuoustime

case 1
    % initialise
    alpha=zeros(model.dim,tseries.ntimes);
    alpha(:,tseries.ntimes)=a(:,tseries.ntimes);

    r=zeros(model.dim,tseries.ntimes);
    N=zeros(model.dim,model.dim,tseries.ntimes);
    L=zeros(model.dim,model.dim,tseries.ntimes);
    V=zeros(model.dim,model.dim,tseries.ntimes);
    V(:,:,tseries.ntimes)=P(:,:,tseries.ntimes);
    % continuous time
    if  strcmp(settings.smoother,'DurbinKoopman')
        
        if settings.multivariate
            
            % backward pass for t=n,...1
            for i=tseries.ntimes:-1:2
                
                if tseries.epochsmissingobs(i)
                    % if there are any missing observations
                    % special treatment for epochs with missing observations:
                    % set Z to 0 (page 111 of Durbin & Koopman)
                    Windex=find(W(i,:));
                    
                    Ztstar=Z(Windex,:,i);% mask out dimensions that are not observed
                    Ftstar=F(Windex,Windex,i);% mask out dimensions that are not observed
                    vtstar=v(Windex,i) ;% mask out dimensions that are not observed
                    
                    L(:,:,i)=T(:,:,i) - K(:,:,i) * Z(:,:,i);%(eq.4)
                    %The weighted sum of innovations 
                    r(:,i-1)=Ztstar' /(Ftstar) * vtstar + L(:,:,i)' * r(:,i);%(eq.2)
                    N(:,:,i-1)=Ztstar' /(Ftstar) * Ztstar + L(:,:,i)' * N(:,:,i) * L(:,:,i);%(eq.3)
                else                 
                    L(:,:,i)=T(:,:,i) - K(:,:,i) * Z(:,:,i);%(eq.4)
                    %The weighted sum of innovations r
                    r(:,i-1)=Z(:,:,i)' /(F(:,:,i)) * v(:,i) + L(:,:,i)' * r(:,i);%(eq.2)
                    N(:,:,i-1)=Z(:,:,i)' /(F(:,:,i)) * Z(:,:,i) + L(:,:,i)' * N(:,:,i) * L(:,:,i);%(eq.3)
                end
                alpha(:,i) = a(:,i) + P(:,:,i) * r(:,i-1);
                V(:,:,i) = P(:,:,i) - P(:,:,i) * N(:,:,i-1) * P(:,:,i);
            end
            
            % separately i=1
            i=1;
            if tseries.epochsmissingobs(i)
                Windex=find(W(i,:));
                Ztstar=Z(Windex,:,i);% mask out dimensions that are not observed
                Ftstar=F(Windex,Windex,i);% mask out dimensions that are not observed
                vtstar=v(Windex,i) ;% mask out dimensions that are not observed
                L(:,:,i)=T(:,:,i) - K(:,:,i) * Z(:,:,i);%(eq.4)
                r0=Ztstar' /(Ftstar) * vtstar + L(:,:,i)' * r(:,i);%(eq.2)
                N0=Ztstar' /(Ftstar) * Ztstar + L(:,:,i)' * N(:,:,i) * L(:,:,i);%(eq.3)
            else
                L(:,:,i)=T(:,:,i) - K(:,:,i) * Z(:,:,i);%(eq.4)
                r0=Z(:,:,i)' /(F(:,:,i)) * v(:,i) + L(:,:,i)' * r(:,i);%(eq.2)
                N0=Z(:,:,i)' /(F(:,:,i)) * Z(:,:,i) + L(:,:,i)' * N(:,:,i) * L(:,:,i);%(eq.3)
            end
            alpha(:,i) = a(:,i) + P(:,:,i) * r0;
            V(:,:,i) = P(:,:,i) - P(:,:,i) * N0 * P(:,:,i);
        else
            % univariate model
            % backward pass for t=n,...1
            
            
            for i=tseries.ntimes:-1:2
                if ~isempty(find(i==tseries.missingepochs))
                    %                special treatment for epochs with missing observations: Z=0
                    %                Koopman & Durbin equation 4.93
                    L(:,:,i)=T(:,:,i);
                    r(:,i-1)=T(:,:,i)'* r(:,i);%(eq.2)
                    N(:,:,i-1)=T(:,:,i)' * N(:,:,i) * T(:,:,i);%(eq.3)    
                else
                    L(:,:,i)=T(:,:,i) - K(:,i) * Z(1,:,i);%(eq.4)
                    r(:,i-1)=Z(1,:,i)' /(F(i)) * v(i) + L(:,:,i)' * r(:,i);%(eq.2)
                    N(:,:,i-1)=Z(1,:,i)' /(F(i)) * Z(1,:,i) + L(:,:,i)' * N(:,:,i) * L(:,:,i);%(eq.3)       
                end
                alpha(:,i) = a(:,i) + P(:,:,i) * r(:,i-1);
                V(:,:,i) = P(:,:,i) - P(:,:,i) * N(:,:,i-1) * P(:,:,i);
            end
            % separately i=1
            i=1;
            L(:,:,i)=T(:,:,i) - K(:,i) * Z(1,:,i);%(eq.4)
            r0=Z(1,:,i)' /(F(i)) * v(i) + L(:,:,i)' * r(:,i);%(eq.2)
            N0=Z(1,:,i)' /(F(i)) * Z(1,:,i) + L(:,:,i)' * N(:,:,i) * L(:,:,i);%(eq.3)
            alpha(:,i) = a(:,i) + P(:,:,i) * r0;
            V(:,:,i) = P(:,:,i) - P(:,:,i) * N0 * P(:,:,i);
        end
        
    elseif strcmp(settings.smoother,'SquareRoot')
    % square root smoother
    % see chapter 6.3 of Durbin & Koopman
    
    if settings.multivariate
        if tseries.missingobs
            error('mssing observations not implemented for square root filter')
        end
        alpha=zeros(model.dim,tseries.ntimes);
        alpha(:,tseries.ntimes)=a(:,tseries.ntimes);

        r=zeros(model.dim,tseries.ntimes);
        N=zeros(model.dim,model.dim,tseries.ntimes);
        V=zeros(model.dim,model.dim,tseries.ntimes);
        V(:,:,tseries.ntimes)=P(:,:,tseries.ntimes);
        Ntilde = zeros(model.dim,model.dim+1,tseries.ntimes);
        % backward pass for t=n,...1
        for i=tseries.ntimes:-1:2
            Ftilde = sqrt(F(:,:,i));
            ustar1 = Ftilde;
            ustar2 = K(:,:,i)*Ftilde;
            L(:,:,i)               = T(:,:,i) - ustar2/(ustar1)*Z(:,:,i);
            Ntmin1star      = [Z(:,:,i)'/(ustar1)' L(:,:,i)'*Ntilde(:,:,i)];
            [Ntplus]        = householder(Ntmin1star);
            Ntilde(:,:,i-1) = Ntplus(1:model.dim,1:model.dim+1);
            N(:,:,i-1)      = Ntilde(:,:,i-1)*Ntilde(:,:,i-1)';
            r(:,i-1)        = Z(:,:,i)' /(F(:,:,i)) * v(:,i) + L(:,:,i)' * r(:,i);%(eq.2)
            alpha(:,i)      = a(:,i) + P(:,:,i) * r(:,i-1);
            V(:,:,i)        = P(:,:,i) - P(:,:,i) * N(:,:,i-1) * P(:,:,i);

        end
        % separately i=1
        i = 1;
        Ftilde = sqrt(F(:,:,i));
        ustar1 = Ftilde;
        ustar2 = K(:,:,i)*Ftilde;
        L(:,:,i)   = T(:,:,i) - ustar2/(ustar1)*Z(:,:,i);
        Ntmin1star = [Z(:,:,i)'*inv(ustar1)' L(:,:,i)'*Ntilde(:,:,i)];
        [Ntplus] = householder(Ntmin1star);
        Ntilde0 = Ntplus(1:model.dim,1:model.dim+1);
        N0         = Ntilde0*Ntilde0';
        r0         = Z(:,:,i)' /(F(:,:,i)) * v(:,i) + L(:,:,i)' * r(:,i);%(eq.2)
        alpha(:,i) = a(:,i) + P(:,:,i) * r0;
        V(:,:,i)   = P(:,:,i) - P(:,:,i) * N0 * P(:,:,i);
    else
        % univariate model
        alpha=zeros(model.dim,tseries.ntimes);
        alpha(:,tseries.ntimes)=a(:,tseries.ntimes);

        r=zeros(model.dim,tseries.ntimes);
        N=zeros(model.dim,model.dim,tseries.ntimes);
        V=zeros(model.dim,model.dim,tseries.ntimes);
        V(:,:,tseries.ntimes)=P(:,:,tseries.ntimes);
        Ntilde = zeros(model.dim,model.dim+1,tseries.ntimes);
        % backward pass for t=n,...1
        for i=tseries.ntimes:-1:2
            Ftilde = sqrt(F(i));
            ustar1 = Ftilde;
            ustar2 = K(:,i)*Ftilde;
            L(:,:,i)               = T(:,:,i) - ustar2/(ustar1)*Z(:,i);
            Ntmin1star      = [Z(:,i)'/(ustar1)' L(:,:,i)'*Ntilde(:,:,i)];
            [Ntplus]        = householder(Ntmin1star);
            Ntilde(:,:,i-1) = Ntplus(1:model.dim,1:model.dim+1);
            N(:,:,i-1)      = Ntilde(:,:,i-1)*Ntilde(:,:,i-1)';
            r(:,i-1)        = Z(:,i)' /(F(i)) * v(i) + L(:,:,i)' * r(:,i);%(eq.2)
            alpha(:,i)      = a(:,i) + P(:,:,i) * r(:,i-1);
            V(:,:,i)        = P(:,:,i) - P(:,:,i) * N(:,:,i-1) * P(:,:,i);

        end
        % separately i=1
        i = 1;
        Ftilde = sqrt(F(i));
        ustar1 = Ftilde;
        ustar2 = K(:,i)*Ftilde;
        L(:,:,i)   = T(:,:,i) - ustar2/(ustar1)*Z(:,i);
        Ntmin1star = [Z(:,i)'*inv(ustar1)' L(:,:,i)'*Ntilde(:,:,i)];
        [Ntplus] = householder(Ntmin1star);
        Ntilde0 = Ntplus(1:model.dim,1:model.dim+1);
        N0         = Ntilde0*Ntilde0';
        r0         = Z(:,i)' /(F(i)) * v(i) + L(:,:,i)' * r(:,i);%(eq.2)
        alpha(:,i) = a(:,i) + P(:,:,i) * r0;
        V(:,:,i)   = P(:,:,i) - P(:,:,i) * N0 * P(:,:,i);
    end
    else
        error('unknown smoother')
    end

case 0
    % integer time steps version
    if  strcmp(settings.smoother,'DurbinKoopman')
        % initialise
        alpha=zeros(model.dim,tseries.ntimes);
        alpha(:,tseries.ntimes)=a(:,tseries.ntimes);

        r=zeros(model.dim,tseries.ntimes);
        N=zeros(model.dim,model.dim,tseries.ntimes);
        L=zeros(model.dim,model.dim,tseries.ntimes);
        V=zeros(model.dim,model.dim,tseries.ntimes);
        V(:,:,tseries.ntimes)=P(:,:,tseries.ntimes);

        % backward pass for t=n,...1

        for i=tseries.ntimes:-1:2
           % 

           if ~isempty(find(i==tseries.missingepochs))

               % special treatment for epochs with missing observations: Z=0
               % Koopman & Durbin equation 4.93
               r(:,i-1)=T'* r(:,i);%(eq.2)
               N(:,:,i-1)=T' * N(:,:,i) * T;%(eq.3)
               alpha(:,i) = a(:,i) + P(:,:,i) * r(:,i-1);
               V(:,:,i) = P(:,:,i) - P(:,:,i) * N(:,:,i-1) * P(:,:,i);
           else
               % normal smoother when observations are available
               L(:,:,i)=T - K(:,i) * Z;%(eq.4)
               r(:,i-1)=Z' /(F(i)) * v(i) + L(:,:,i)' * r(:,i);%(eq.2)
               N(:,:,i-1)=Z' /(F(i)) * Z + L(:,:,i)' * N(:,:,i) * L(:,:,i);%(eq.3)
               alpha(:,i) = a(:,i) + P(:,:,i) * r(:,i-1);
               V(:,:,i) = P(:,:,i) - P(:,:,i) * N(:,:,i-1) * P(:,:,i);
           end
        end
        % separately i=1
        i=1;
        L(:,:,i)=T - K(:,i) * Z;%(eq.4)
        r0=Z' /(F(i)) * v(i) + L(:,:,i)' * r(:,i);%(eq.2)
        N0=Z' /(F(i)) * Z + L(:,:,i)' * N(:,:,i) * L(:,:,i);%(eq.3)
        alpha(:,i) = a(:,i) + P(:,:,i) * r0;
        V(:,:,i) = P(:,:,i) - P(:,:,i) * N0 * P(:,:,i);

    else
        error('unknown smoother')
    end

end

% save in structure
estsmooth.alpha=alpha;
estsmooth.V=V;
estsmooth.L=L;
estsmooth.N=N;

end