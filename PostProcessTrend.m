function [meanslope] = PostProcessTrend(settings,tseries,est,estvar,estcov,model)
% PostProcessTrend.m computes a mean linear trend and its error variance
% based on the stochastic trend and its error variance.
%
% HOW:     [meanslope] = PostProcessTrend(tseries,est,estvar,covsmooth)
%
% Input:
%       settings, of which this function uses:
%           normalised
%
%       tseries structure which should contain at least:
%           time [epochs x 1]           unnormalized time
%           ntimes                      number of epochs
%
%       est, structure of which this function uses:
%           est,slope [epochs x 1]      time variable slope
%
%       estvar, structure of which this function uses:
%           estvar.slope [epochs x 1]   slope error variances%
%
%       covsmooth, structure of which this function uses:
%           estcov.alpha                matrix cov(dim,dim,cov_alpha_t,cov_alpha_j)
%              autocovariance of (smoothed) states
%
%           model                           contains indices of components
%
% Output:
%       meanslope, structure containing:
%           meanslope.slope                   mean slope
%           meanslope.slopevar                error variance of mean slope
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
% Version 1.0 August 2014 DBT Broerse
%
% Version 1.1 March 2015 DBT Broerse
% - inclusion of correlations between levels and slopes at different epochs
%
% Version 1.2 April 2015 DBT Broerse
% - option for non normalised time series
%
% Version 1.3 June 2018 DBT Broerse
% - extension to multivariate models
%
% Version 1.4 April 2020 DBT Broerse
% - fixed bug mean slope variance
% - divided time step by two, instead of correction for twice time step at
% the end
%
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
%----------------------------------------------------------------------------

narginchk(6,6);

meanslope.slope=[];
meanslope.slopevar=[];


if settings.slope
% total time is time between halfway last two points minus halfway time of
% first two points
totaltime=tseries.time(end)-tseries.time(1);
% mean slope is mean slope per time step integrated over time, divided by
% total time

% dt is time step until half way the interval
% this definition is consistent with the total time
dt(1)=(tseries.time(2)-tseries.time(1))*0.5;
for i=2:tseries.ntimes-1
    dt(i)=0.5*(tseries.time(i+1)-tseries.time(i-1));
end
% last
dt(tseries.ntimes)=(tseries.time(tseries.ntimes)-tseries.time(tseries.ntimes-1))*0.5;



for p=1:tseries.ntseries
    % the mean of the slope
    % could also be estimated based on first and last step, see Frederikse
    % et al. 2016, but for the uncertainty propagation, the mean has to be
    % taken properly
    meanslope.slope{p}=1/(totaltime)*sum(est.slope{p}.*dt);

    % variance
    var_contr=sum(estvar.slope{p}.*dt.^2);
    
    % now add covariances
    % I'm sure this can be optimized
    cov_contr=0;
    if settings.normalise
        for i=1:tseries.ntimes
            for j=1:tseries.ntimes
                if j~=i
                    normconstant=1/tseries.period^2;
                    % estcov has normalised values
                    cov_contr=cov_contr+dt(i)*dt(j)*estcov.alpha(model.indexslope{p},model.indexslope{p},i,j)*normconstant;
                end
            end
        end
    else
        for i=1:tseries.ntimes
            for j=1:tseries.ntimes
                if j~=i
                    cov_contr=cov_contr+dt(i)*dt(j)*estcov.alpha(model.indexslope{p},model.indexslope{p},i,j);
                end
            end
        end
    end
    
    % add variance and covariance contributions
    meanslope.slopevar{p}=1/(totaltime^2)*(var_contr+cov_contr);
    
end

% %% approach from Frederikse et al. 2015 JGR 
% % based on begin and end trend
% meanslope.slope_begin_end=[];
% meanslope.slopevar_begin_end=[];
% for p=1:tseries.ntseries
%     % begin and end trend, equation 43 Frederikse et al.
%     meanslope.slope_begin_end{p} = 1/totaltime * ( est.trend{p}(end) - est.trend{p}(1));
%     
%     % variance
%    % no normalisation is needed, as it concerns the level, not slope, whihc
%     % is time independent
%     covar = estcov.alpha(model.indexlevel{p},model.indexlevel{p},tseries.ntimesobs,1);
%  
%     % see equation 44 Frederikse et al.
%     meanslope.slopevar_begin_end{p} = 1/(totaltime^2) * (estvar.slope{p}(end) + estvar.slope{p}(1) - 2*covar);
% 
% end
end

end

