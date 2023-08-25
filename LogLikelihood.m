function [LogL]= LogLikelihood(tseries,estfilter,settings,model)

% LogLikelihood.m calculates the log of the likelihood L(Y_n).
% Here L(Y_n) = p(y_1,...,y_n) = p(y_1) PRODUCT_t=2:n ( p(y_t|Y_t-1) )
%
%
% HOW:     [LogL]= LogLikelihood(tseries,estfilter,settings)
%
% Input:    
%       tseries, structure of which this function uses:
%           tseries.t [epochs x 1]          normalised time vector
%           tseries.ntimes                  number of epochs
%           tseries.missingepochs           indexes of missing data (only for integer time)
%           tseries.np                      total number of single observations
%
%       estsmooth, structure that contains
%           estfilter.v [epochs]            Prediction error
%           estfilter.F [epochs]            The variance of prediction error v
%   
% Output:   LogL                            log likelihood
%
%
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
% uses: none
% 
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 June 2014 DBT Broerse
%
% Version 1.1 August 2014 DBT Broerse
% - correction for log(abs(F)) when observation is missing
% - correction for ntimes when epochs are missing
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
%
% Version 2.1 June 2018 DBT Broerse
% - extended to multivariate models
%
% Version 2.2 May 2019 DBT Broerse
% - correction of log likelihood equation
%  % previous version contained the absolute instead of
% determinant of F 
%
% Version 2.3 June 2019 DBT Broerse
% - missing observations in multivariate models fixes
%
%----------------------------------------------------------------------------
% remarks: 
%----------------------------------------------------------------------------
%
%%%% THEORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Formulation of the log likelihood according to chapter 7 of Durbin and
% Koopman, equation 7.2.
%
% log L (Y_n) = -ntimes/2 log(2pi) - 1/2 sum_t=1:n (log(det(F_t)) +
% v_t'*inv(F_t)*v_t)
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract from structure

F=estfilter.F;
v=estfilter.v;

if settings.multivariate
    W=model.W;
end

LogL=[];

%% compute log likelihood

% take sum of v_t'*inv(F_t)*v_t)


switch settings.multivariate
    
    case 0
        
        % univariate case
        SumTemp=0;SumTemp2=0;
        for i=1:tseries.ntimes
            % skip epochs without observations
            if isempty(find(i==tseries.missingepochs))
               % SumTemp=SumTemp+v(i)'*inv(F(i))*v(i);
                SumTemp=SumTemp+v(i)'*(F(i)\v(i));
                SumTemp2=SumTemp2+log(det(F(i)));
%             else
%                 i
%                 disp('skip')
            end
            
        end
        
   
    case 1
        
        % multivariate case
        SumTemp=0;SumTemp2=0;
        for i=1:tseries.ntimes

           % check for missing epochs, because v and F are defined as zero for
           % missing epochs
               % SumTemptemp=SumTemp+v(:,i)'*inv(F(:,:,i))*v(:,i);
               if tseries.epochsmissingobs(i)
                   % only take observed dimensions
                   Windex=find(W(i,:));
                   Ftstar=F(Windex,Windex,i);% mask out dimensions that are not observed
                   vtstar=v(Windex,i) ;% mask out dimensions that are not observed
                   SumTemp=SumTemp+vtstar'*(Ftstar\vtstar);
               else
                   vt=v(:,i);
                   Ft=F(:,:,i);
                   SumTemp=SumTemp+vt'*(Ft\vt);
               end
                
               
                % previous version contained the absolute instead of
                % determinant
                if det(F(:,:,i)) < 0
                    warning('correlation > 1 in F (one-step-ahead prediction error variance)')
                    % reset to zero
                   F(:,:,i)=zeros(size(F(:,:,i))); 
                end
                
                if tseries.epochsmissingobs(i)
                    SumTemp2=SumTemp2+log(det(Ftstar));
                else
                    SumTemp2=SumTemp2+log(det(Ft));
                end
                
            
        end
       
end



% log likelihood
LogL = -1/2*(tseries.np*log(2*pi)+SumTemp2+SumTemp);

if ~isreal(LogL)
    SumTemp2
    F(:,:,end)
    det(F(:,:,end))
    log(det(F(:,:,end)))
    SumTemp
    disp('log likelihood is complex')
    LogL=-1e3;
end


end




