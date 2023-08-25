function [estsmooth] = DisturbanceSmoother(tseries,model,settings,estfilter,estsmooth)

% Disturbance smoother computes the smoothed disturbance vectors epsilon
% (disturbance) and eta (process noise)
%
%
% HOW:     [estsmooth] = DisturbanceSmoother(tseries,model,settings,estfilter,estsmooth)
%
% Input:
%       tseries, structure of which this function uses:
%           tseries.t [epochs x 1]      normalised time vector
%           tseries.ntimes              length of normalised time vector
%           tseries.missingepochs       indexes with missing data (integer time only)
%
%       model, structure of which this function uses:
%           model.Z   [epochs x n]  design matrix
%           model.T   [n x n x epochs]  transition matrix
%           model.Q   [n x n x epochs]  process noise matrix
%           model.H   [epochs x epochs] observational noise matrix
%           model.R   [n x n]           selection matrix
%
%       estfilter, structure of which this function uses:
%           estfilter.v   [epochs x 1]  prediction error
%           estfilter.K   [n x epochs]            Kalman gain
%           estfilter.F   [epochs x epochs]       The variance of prediction error v_t
%
%       settings, structure of which this function uses:
%           settings.distsmoother       type of smoother
%                                       - 'DurbinKoopman' according to Durbin & Koopman
%
%       estsmooth (only when estsmooth.alpha and estsmooth.V already exist)
%
% Output:
%
%       estsmooth, structure that contains
%           estsmooth.epsilon           smoothed disturbance
%           estsmooth.eta               smoothed process noise
%           estsmooth.varepsilon        variance of disturbance given Yn
%           estsmooth.vareta=Vareta     variance of process noise given Yn
%           estsmooth.u                 smoothing error
%           estsmooth.D                 smoothing error variance
%           estsmooth.r                 weighted sum of innovations
%           estsmooth.N                 variance of weighted sum of innovations
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
% Version 2.1 March 2015 DBT Broerse
% - changed size of Z matrix
%
% Version 2.2 June 2018 DBT Broerse
% - extended to multivariate models
% - changed indexing of Z for continuous time
%
% Version 2.3 June 2019 DBT Broerse
% - addition of missing epochs in multivariate models
%
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------
%
%
%%%% THEORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Formulation of the disturbance smoothing according to chapter 4.5.3 of Durbin and
% Koopman, equation 4.69.
%
% Disturbance smoothing consist of following recursion for t=n,..,1
%
% smoothing error:
%
% u_t = inv(F_t) v_t - K_t'r_t                                       (eq.1)
%
% D_t = inv(F_t) + K_t' N_t K_t                                      (eq.2)
%
% weighted sum of the inverse variances of innovations N, with N_n = 0
%
% N_t-1 = Z_t' D_t Z_t + T_t' N_t T_t - Z_t' K_t' N_t T_t - T_t' N_t K_t Z_t (eq.3)
%
% weighted sum of innovations, with r_n = 0
%
% r_t-1 = Z_t' u_t + T_t' r_t                                        (eq.4)
%
% smoothed disturbance and its conditional variance:
%
% epsilon_t = H_t u_t                                                (eq.5)
% Var_epsilon_t = H_t - H_t D_t H_t
%
% smoothed process noise and its conditional variance:
%
% eta_t = Q_t R_t' r_t
% Var_eta_t = Q_t - Q_t R_t' N_t' R_t Q_t                            (eq.6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
narginchk(4,5);
%error(nargchk(4,5,nargin));

% set default smoother type
if ~any(strcmp('distsmoother',fieldnames(settings)))
    settings.distsmoother=[];
end
if isempty(settings.distsmoother)
    if isempty(settings.smoother)
        if isempty(settings.filter)
            settings.distsmoother='DurbinKoopman';
        else
            settings.distsmoother=settings.filter;
        end
    else
        settings.distsmoother=settings.smoother;
    end
end

% extract form structure
K=estfilter.K;
F=estfilter.F;
v=estfilter.v;

Z=model.Z;
Q=model.Q;
R=model.R;
T=model.T;
H=model.H;
if settings.multivariate
W=model.W;
end

% different code for continuous time (general case) and integer time steps
switch settings.continuoustime
    
    case 1
        
        
        
        
        % initialise
        epsilon=zeros(tseries.ntseries,tseries.ntimes);
        eta=zeros(model.dim,tseries.ntimes);
        Vareta=zeros(model.dim,model.dim,tseries.ntimes);
        Varepsilon=zeros(tseries.ntseries,tseries.ntseries,tseries.ntimes);
        u=zeros(tseries.ntseries,tseries.ntimes);
        D=zeros(tseries.ntseries,tseries.ntseries,tseries.ntimes);
        r=zeros(model.dim,tseries.ntimes);
        N=zeros(model.dim,model.dim,tseries.ntimes);
        Ntilde = zeros(model.dim,model.dim+1,tseries.ntimes);
        
        % remove redundant dimensions for univariate models
        if ~settings.multivariate
            Varepsilon=squeeze(Varepsilon);
            D=squeeze(D);
        end
        
        if  strcmp(settings.distsmoother,'DurbinKoopman')
            % backward pass for t=n,...1
            if settings.multivariate
                for i=tseries.ntimes:-1:2
                    if tseries.epochsmissingobs(i)
                        % if there are any missing observations
                        % special treatment for epochs with missing observations:
                        % set Z to 0 (page 111 of Durbin & Koopman)
                        Windex=find(W(i,:));
                        Ftstar=F(Windex,Windex,i);% mask out dimensions that are not observed
                        vtstar=v(Windex,i) ;% mask out dimensions that are not observed
                        
                        u(Windex,i)=(Ftstar)\vtstar-K(:,Windex,i)'*r(:,i); %(eq.1)
                        D(Windex,Windex,i)=inv(Ftstar)+K(:,Windex,i)'*N(:,:,i)*K(:,Windex,i); %(eq.2) 
                    else
                        u(:,i)=(F(:,:,i))\v(:,i)-K(:,:,i)'*r(:,i); %(eq.1)
                        D(:,:,i)=inv(F(:,:,i))+K(:,:,i)'*N(:,:,i)*K(:,:,i); %(eq.2)
                    end
                    
                    N(:,:,i-1)=Z(:,:,i)'*D(:,:,i)*Z(:,:,i) + T(:,:,i)'*N(:,:,i)*T(:,:,i) -Z(:,:,i)'*K(:,:,i)'*N(:,:,i)*T(:,:,i)-T(:,:,i)'*N(:,:,i)*K(:,:,i)*Z(:,:,i);%(eq.3)
                    r(:,i-1)=Z(:,:,i)'*u(:,i)+T(:,:,i)'*r(:,i);%(eq.4)
                    epsilon(:,i)=H(:,:,i)*u(:,i);%(eq.5)
                    Varepsilon(:,:,i)=H(:,:,i)-H(:,:,i)*D(:,:,i)*H(:,:,i);%(eq.5)
                    eta(:,i)=Q(:,:,i)*R'*r(:,i);%(eq.6)
                    Vareta(:,:,i)=Q(:,:,i)-Q(:,:,i)*R'*N(:,:,i)*R*Q(:,:,i);%(eq.6)
                    
                    
                end
                i=1;
                if tseries.epochsmissingobs(i)
                    % if there are any missing observations
                    % special treatment for epochs with missing observations:
                    % set Z to 0 (page 111 of Durbin & Koopman)
                    Windex=find(W(i,:));
                    Ftstar=F(Windex,Windex,i);% mask out dimensions that are not observed
                    vtstar=v(Windex,i) ;% mask out dimensions that are not observed
                    
                    u(Windex,i)=(Ftstar)\vtstar-K(:,Windex,i)'*r(:,i); %(eq.1)
                    D(Windex,Windex,i)=inv(Ftstar)+K(:,Windex,i)'*N(:,:,i)*K(:,Windex,i); %(eq.2)
                else
                    u(:,i)=(F(:,:,i))\v(:,i)-K(:,:,i)'*r(:,i); %(eq.1)
                    D(:,:,i)=inv(F(:,:,i))+K(:,:,i)'*N(:,:,i)*K(:,:,i); %(eq.2)
                end
                epsilon(:,i)=H(:,:,i)*u(:,i);%(eq.5)
                Varepsilon(:,:,i)=H(:,:,i)-H(:,:,i)*D(:,:,i)*H(:,:,i);%(eq.5)
                eta(:,i)=Q(:,:,i)*R'*r(:,i);%(eq.6)
                Vareta(:,:,i)=Q(:,:,i)-Q(:,:,i)*R'*N(:,:,i)*R*Q(:,:,i);%(eq.6)
                
            else
                
                % univariate version
                for i=tseries.ntimes:-1:2
                    % skip epochs without observations
                    if isempty(find(i==tseries.missingepochs,1))
                        u(i)=(F(i))\v(i)-K(:,i)'*r(:,i); %(eq.1)
                        D(i)=(F(i))\1+K(:,i)'*N(:,:,i)*K(:,i); %(eq.2)
                        N(:,:,i-1)=Z(1,:,i)'*D(i)*Z(1,:,i) + T(:,:,i)'*N(:,:,i)*T(:,:,i) -Z(1,:,i)'*K(:,i)'*N(:,:,i)*T(:,:,i)-T(:,:,i)'*N(:,:,i)*K(:,i)*Z(1,:,i);%(eq.3)
                        r(:,i-1)=Z(1,:,i)'*u(i)+T(:,:,i)'*r(:,i);%(eq.4)
                        epsilon(i)=H(i,i)*u(i);%(eq.5)
                        Varepsilon(i)=H(i,i)-H(i,i)*D(i)*H(i,i);%(eq.5)
                        eta(:,i)=Q(:,:,i)*R'*r(:,i);%(eq.6)
                        Vareta(:,:,i)=Q(:,:,i)-Q(:,:,i)*R'*N(:,:,i)*R*Q(:,:,i);%(eq.6)
                    end
                    
                end
                i=1;
                u(i)=(F(i))\v(i)-K(:,i)'*r(:,i); %(eq.1)
                D(i)=(F(i))\1+K(:,i)'*N(:,:,i)*K(:,i); %(eq.2)
                epsilon(i)=H(i,i)*u(i);%(eq.5)
                Varepsilon(i)=H(i,i)-H(i,i)*D(i)*H(i,i);%(eq.5)
                eta(:,i)=Q(:,:,i)*R'*r(:,i);%(eq.6)
                Vareta(:,:,i)=Q(:,:,i)-Q(:,:,i)*R'*N(:,:,i)*R*Q(:,:,i);%(eq.6)
            end
        elseif strcmp(settings.distsmoother,'SquareRoot')
            % see chapter 6.3 of Durbin & Koopman
            
            
            if settings.multivariate
                for i=tseries.ntimes:-1:2
                    u(:,i)=(F(:,:,i))\v(:,i)-K(:,:,i)'*r(:,i); %(eq.1)
                    D(:,:,i)=inv(F(:,:,i))+K(:,:,i)'*N(:,:,i)*K(:,:,i); %(eq.2)
                    ustar1 = sqrt(F(:,:,i));
                    ustar2 = K(:,:,i)*ustar1;
                    L          = T(:,:,i) - ustar2/(ustar1)*Z(:,:,i);
                    Ntmin1star = [Z(:,:,i)'*inv(ustar1)' L'*Ntilde(:,:,i)];
                    [Ntplus] = householder(Ntmin1star);
                    Ntilde(:,:,i-1) = Ntplus(1:model.dim,1:model.dim+1);
                    N(:,:,i-1) = Ntilde(:,:,i-1)*Ntilde(:,:,i-1)';
                    r(:,i-1)=Z(:,:,i)'*u(:,i)+T(:,:,i)'*r(:,i);%(eq.4)
                    epsilon(:,i)=H(:,:,i)*u(:,i);%(eq.5)
                    Varepsilon(:,:,i)=H(:,:,i)-H(:,:,i)*D(:,:,i)*H(:,:,i);%(eq.5)
                    eta(:,i)=Q(:,:,i)*R'*r(:,i);%(eq.6)
                    Vareta(:,:,i)=Q(:,:,i)-Q(:,:,i)*R'*N(:,:,i)*R*Q(:,:,i);%(eq.6)
                end
                i=1;
                u(:,i)=(F(:,:,i))\v(:,i)-K(:,:,i)'*r(:,i); %(eq.1)
                D(:,:,i)=inv(F(:,:,i))+K(:,:,i)'*N(:,:,i)*K(:,:,i); %(eq.2)
                epsilon(:,i)=H(:,:,i)*u(:,i);%(eq.5)
                Varepsilon(:,:,i)=H(:,:,i)-H(:,:,i)*D(:,:,i)*H(:,:,i);%(eq.5)
                eta(:,i)=Q(:,:,i)*R'*r(:,i);%(eq.6)
                Vareta(:,:,i)=Q(:,:,i)-Q(:,:,i)*R'*N(:,:,i)*R*Q(:,:,i);%(eq.6)
            else
                % univariate model
                for i=tseries.ntimes:-1:2
                    u(i)=(F(i))\v(i)-K(:,i)'*r(:,i); %(eq.1)
                    D(i)=inv(F(i))+K(:,i)'*N(:,:,i)*K(:,i); %(eq.2)
                    ustar1 = sqrt(F(i));
                    ustar2 = K(:,i)*ustar1;
                    L          = T(:,:,i) - ustar2/(ustar1)*Z(1,:,i);
                    Ntmin1star = [Z(1,:,i)'*inv(ustar1)' L'*Ntilde(:,:,i)];
                    [Ntplus] = householder(Ntmin1star);
                    Ntilde(:,:,i-1) = Ntplus(1:model.dim,1:model.dim+1);
                    N(:,:,i-1) = Ntilde(:,:,i-1)*Ntilde(:,:,i-1)';
                    r(:,i-1)=Z(1,:,i)'*u(i)+T(:,:,i)'*r(:,i);%(eq.4)
                    epsilon(i)=H(i,i)*u(i);%(eq.5)
                    Varepsilon(i)=H(i,i)-H(i,i)*D(i)*H(i,i);%(eq.5)
                    eta(:,i)=Q(:,:,i)*R'*r(:,i);%(eq.6)
                    Vareta(:,:,i)=Q(:,:,i)-Q(:,:,i)*R'*N(:,:,i)*R*Q(:,:,i);%(eq.6)
                end
                i=1;
                u(i)=(F(i))\v(i)-K(:,i)'*r(:,i); %(eq.1)
                D(i)=inv(F(i))+K(:,i)'*N(:,:,i)*K(:,i); %(eq.2)
                epsilon(i)=H(i,i)*u(i);%(eq.5)
                Varepsilon(i)=H(i,i)-H(i,i)*D(i)*H(i,i);%(eq.5)
                eta(:,i)=Q(:,:,i)*R'*r(:,i);%(eq.6)
                Vareta(:,:,i)=Q(:,:,i)-Q(:,:,i)*R'*N(:,:,i)*R*Q(:,:,i);%(eq.6)
            end
        else
            error('unknown smoother')
        end
        
        
    case 0
        
        % integer time steps
        if  strcmp(settings.distsmoother,'DurbinKoopman')
            % initialise
            epsilon=zeros(1,tseries.ntimes);
            eta=zeros(model.dim,tseries.ntimes);
            Vareta=zeros(model.dim,model.dim,tseries.ntimes);
            Varepsilon=zeros(tseries.ntimes,1);
            u=zeros(tseries.ntimes,1);
            D=zeros(tseries.ntimes,1);
            
            r=zeros(model.dim,tseries.ntimes);
            N=zeros(model.dim,model.dim,tseries.ntimes);
            
            % backward pass for t=n,...1
            
            for i=tseries.ntimes:-1:2
                % skip epochs without observations
                if isempty(find(i==tseries.missingepochs))
                    u(i)=(F(i))\v(i)-K(:,i)'*r(:,i); %(eq.1)
                    D(i)=inv(F(i))+K(:,i)'*N(:,:,i)*K(:,i); %(eq.2)
                    N(:,:,i-1)=Z'*D(i)*Z + T'*N(:,:,i)*T -Z'*K(:,i)'*N(:,:,i)*T-T'*N(:,:,i)*K(:,i)*Z;%(eq.3)
                    r(:,i-1)=Z'*u(i)+T'*r(:,i);%(eq.4)
                    epsilon(1,i)=H(i,i)*u(i);%(eq.5)
                    Varepsilon(i)=H(i,i)-H(i,i)*D(i)*H(i,i);%(eq.5)
                    eta(:,i)=Q*R'*r(:,i);%(eq.6)
                    Vareta(:,:,i)=Q-Q*R'*N(:,:,i)*R*Q;%(eq.6)
                end
                
            end
            i=1;
            u(i)=(F(i))\v(i)-K(:,i)'*r(:,i); %(eq.1)
            D(i)=inv(F(i))+K(:,i)'*N(:,:,i)*K(:,i); %(eq.2)
            epsilon(1,i)=H(i,i)*u(i);%(eq.5)
            Varepsilon(i)=H(i,i)-H(i,i)*D(i)*H(i,i);%(eq.5)
            eta(:,i)=Q*R'*r(:,i);%(eq.6)\
            Vareta(:,:,i)=Q-Q*R'*N(:,:,i)*R*Q;%(eq.6)
        else
            error('unknown smoother')
        end
        
end

% save in structure
estsmooth.epsilon=epsilon;
estsmooth.eta=eta;
estsmooth.varepsilon=Varepsilon;
estsmooth.vareta=Vareta;
estsmooth.u=u;
estsmooth.D=D;
estsmooth.r=r;
estsmooth.N=N;

end

