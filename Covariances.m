function [covsmooth] = Covariances(tseries,estfilter,estsmooth,model)
% Covariances computes the covariance matrices for smoothed estimates.
%
% HOW:     [covsmooth] = Covariances(tseries,estfilter,estsmooth,model)
%
% Input:
%       tseries structure which should contain at least:
%           ntimes                      number of epochs
%
%           model                           contains indices of components
%
%
% Output:   
%       cov, structure containing:
%           covsmooth.alpha                    matrix cov(dim,dim,alpha_t,alpha_j)
%           covsmooth.corr_coef                matrix cov(dim,dim,corr_alpha_t,corr_alpha_j)
%              provides Pearson correlation coefficient                            
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
% Covariance is for alpha_i,alpha_j is computed as:
%
% Cov(alpha_t,alpha_j) = P_t * L'_t * L'_t+1 .. * L'_j-1 * ( I - N_j-1*P_j)
%
% see table 4.4 of Durbin and Koopman (2012)
%
% correlation coefficient is defined as:
%
% corrcoef_t,j = cov_t,j / ( stdev_t * stdev_j)
%
%----------------------------------------------------------------------------
% uses: none
% 
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 October 2014 DBT Broerse
%
% Version 1.1 March 2015 DBT Broerse
% Fixed a few errors in the computation of cov(alpha,alpha)
%
%
%----------------------------------------------------------------------------
% remarks: currently only cov(alpha_t,alpha_j) is computed.
% loops can be set up more efficient
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
%----------------------------------------------------------------------------

% initialise
covsmooth.alpha=zeros(model.dim,model.dim,tseries.ntimes,tseries.ntimes);

P=estfilter.P;
L=estsmooth.L;
N=estsmooth.N;
V=estsmooth.V;
% identity matrix
I=eye(model.dim,model.dim);

L_transp=zeros(size(L));
for t=1:tseries.ntimes
    L_transp(:,:,t)=L(:,:,t)';
end
%return
% covariance of [alpha_i,alpha_j] for j >= t
for t=1:tseries.ntimes
    % L'_{t+1} .. L'_{j-1} becomes Identity matrix for j=t+1
    if t<tseries.ntimes
        j=t+1;
        covsmooth.alpha(:,:,t,j)=P(:,:,t)*L_transp(:,:,t)*(I-N(:,:,j-1)*P(:,:,j));
    end
    for j=t+2:tseries.ntimes
         % initialise the product L'_t * L'_t+1 .. * L'_j-1
         L_prod=I; 
         for k=t:j-1
              L_prod=L_prod*L_transp(:,:,k);
         end
         covsmooth.alpha(:,:,t,j)=P(:,:,t)*L_prod*(I-N(:,:,j-1)*P(:,:,j));
       
   end    
end

% for t=1:tseries.ntimes
%    for j=1:tseries.ntimes
%        if (j==t+1)
%            covsmooth.alpha(:,:,t,j)=P(:,:,t)*L_transp(:,:,t)*(I-N(:,:,j-1)*P(:,:,j));
%            %return
%        elseif (j>t+1)   
%          % initialise the product L'_t * L'_t+1 .. * L'_j-1
%          L_prod=I; 
%          for k=t:j-1
%               L_prod=L_prod*L_transp(:,:,k);
%          end
%          covsmooth.alpha(:,:,t,j)=P(:,:,t)*L_prod*(I-N(:,:,j-1)*P(:,:,j));
%        end
%    end    
% end

% now copy contents to alpha for t<j 

for t=1:tseries.ntimes
   for j=1:tseries.ntimes
       if (j < t)
           % covariance matrices are not diagonal, so use transpose
           % as cov(alpha_i,alpha_j) == cov(alpha_j,alpha_t)^T
            covsmooth.alpha(:,:,t,j)=covsmooth.alpha(:,:,j,t)';
       end
   end    
end


%% compute correlation coefficient
covsmooth.corr_coef=zeros(size(covsmooth.alpha));

% first compute standard deviations
stdev=zeros(model.dim,tseries.ntimes);
for i=1:tseries.ntimes
    stdev(:,i)=sqrt(diag(V(:,:,i)));
end

% likely this can be done without loop
for t=1:tseries.ntimes
   for j=1:tseries.ntimes
       % check if this is correct
      covsmooth.corr_coef(:,:,t,j)=covsmooth.alpha(:,:,t,j)./(stdev(:,t)*stdev(:,j)');
  
%    % next lines leads to same results for correlations
%       for a=1:model.dim
%           for b=1:model.dim
%             covsmooth.corr_coef_check(a,b,t,j)=covsmooth.alpha(a,b,t,j)/...
%                 sqrt(V(a,a,t)*V(b,b,j));
%           end
%       end
   end
end

end

