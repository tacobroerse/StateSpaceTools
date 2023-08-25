function [tseries,settings]=StateSpacePreProcess(tseries,settings)
% StateSpacePreProcess normalises time (leading to t=1,2,..n when equal time
% steps). 
%
% HOW:     [tseries,settings]=StateSpacePreProcess(tseries,settings)
%
% Input:    tseries structure which should contain at least:
%               tseries.time [epochs x 1]           time vector
%               tseries.y [epochs x 1]              observation vector 
%               tseries.missingepochs               indexes of 'time' with missing data (for integer time only) 
%           
%           settings structure which should contain at least:
%           	settings.checktime                   whether program should check whether time is indeed equally spaced
%           	settings.continuoustime              user provided type of time steps,
%                                                       if equally spaced, continuoustime is false
%               settings.periodscycle                periods of cylce terms ( [] if there are no cycle components)
%           	settings.interventiontime            time of intervention (step function)
%               settings.normalise                   whether time should be
%               normalized (only possible if continuoustime == true)
%
%
% Output:   
%           tseries, with added contents
%               tseries.t                           normalised time
%               tseries.Y                           observation vector with NaN's at
%                                                   missing epochs
%               tseries.period                      observational period
%               tseries.frequency                   measurement frequency (average for
%                                                   continuous time)
%               tseries.datagaps                    true if there is missing data and
%                                                   time steps are equally spaced
%               tseries.indexavailable              index of available observations
%
%           settings, with added contents
%               settings.numbercycles               number of cycle terms
%               settings.lambda                     normalised cycle frequency
%               settings.tau                        time index, startpoint for step function
%               settings.numberinterventions        number of interventions
%               settings.numberregressors           number of regressors           
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
% - changes in comments
%
% Version 1.2 September 2014 DBT Broerse
% - added regression
% - added step intervention
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
%
% Version 2.1 April 2015 DBT Broerse
% - normalisation is now optional
%
% Version 2.2 May 2015 DBT Broerse
% - shift to left of tseries.dt
%
% Version 2.3 June 2018 DBT Broerse
% - extended to multivariate models
%
%
%----------------------------------------------------------------------------
% remarks: 
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
%----------------------------------------------------------------------------
narginchk(2,2);

% extract from structure

% index of available observations 

tseries.datagaps=false;

% do some checks on the inputs
[rows cols] =size(tseries.time);
if (rows == 1) 
    tseries.time=tseries.time';
   %error('time should be organised in rows')
end

tseries.ntimes=length(tseries.time);

% check ordering of y
[RowsY,ColsY]=size(tseries.y);
if RowsY == tseries.ntimes
    tseries.y=tseries.y';
    tseries.ntseries = ColsY;
elseif ColsY == tseries.ntimes
    tseries.ntseries = RowsY;
else
    disp('error, observation vector has different length than time vector')
    tseries.ntimes
    nobs
    return
end
    
if settings.AR
    if settings.ARorder > 1
        disp('strijk je gaffel mop')
        error('Q in StateSpaceSetup and ExtractVarainceDisturbance are not yet correct for AR orders > 1')
    end
end




if ~settings.multivariate && tseries.ntseries > 1
    error('settings.multivariate set to false while tseries.ntseries is larger than one')
elseif tseries.ntseries == 1 && settings.multivariate
    disp('found only one time series, resetting to univariate model')
    settings.multivariate = 0;
    tseries.missingepochs = tseries.missingepochs{1};
end

    

settings.lambda=[];settings.tau=[];
if settings.normalise
    if strcmp(settings.normalisetimeby,'mean')
    % (Average) time between measurements
    tseries.period=(tseries.time(end)-tseries.time(1))/(tseries.ntimes-1);
    elseif strcmp(settings.normalisetimeby,'median')
        tseries.period=median(diff(tseries.time));
    else
        error('time normalisation not set properly')
    end
    % Sampling frequency (average in case of inequally space time)
    tseries.frequency=1/tseries.period;

    % Normalized time (integers in case of equally spaced time)
    tseries.t=(tseries.time-tseries.time(1))/tseries.period+1;

    % Normalized intervention time and index
    if ~isempty(settings.interventiontime)
        settings.numberinterventions=length(settings.interventiontime);
        for i=1:settings.numberinterventions
            % check whether the intervention time lies within time series
            % time
            
            InterventionTimeNorm(i)=(settings.interventiontime(i)-tseries.time(1))/tseries.period+1;
            if InterventionTimeNorm(i) < tseries.t(end)
            
            % check which times are above this number
            settings.tau(i)=find(tseries.t>InterventionTimeNorm(i),1,'first');
            else
                disp(strcat('WARNING: intervention time:',num2str(settings.interventiontime(i)),' exceeds available time'))
                disp('reducting number of interventions')
                settings.numberinterventions=settings.numberinterventions-1;
                settings.interventiontime(i)=[];
            end
        end
    else
        settings.numberinterventions=0;
    end
    
   
else
    if ~settings.continuoustime
        error('option normalize false only possible for continuoustime')
    end
    % (Average) time between measurements
    tseries.period=(tseries.time(end)-tseries.time(1))/(tseries.ntimes-1);
    % Sampling frequency (average in case of inequally space time)
    tseries.frequency=1/tseries.period;

    % Normalized time (integers in case of equally spaced time)
    tseries.t=tseries.time;

    % Normalized intervention time and index
    if ~isempty(settings.interventiontime)
        settings.numberinterventions=length(settings.interventiontime);
        for i=1:settings.numberinterventions
            InterventionTimeNorm(i)=settings.interventiontime;
            % check which times are above this number
            settings.tau(i)=find(tNorm>InterventionTimeNorm(i),1,'first');
        end
    else
        settings.numberinterventions=0;
    end    
end



% last time step is median time step
% time step should be one ahead, since it is used for one-step-ahead
% prediction
tseries.dt(1:tseries.ntimes-1)=diff(tseries.t);
% average time step (differs from tseries.period in the sense that the
% latter does not take normalised time into account)
tseries.median_dt=median(tseries.dt(1:tseries.ntimes-1));
% set tseries.dt(1) to the median time step
tseries.dt(tseries.ntimes)=tseries.median_dt;


% number of regressors
if (settings.regression) 
    [settings.numberregressors dummy]=size(tseries.regressors);
    % check whether size of regressor is right, or that it has to be
    % transposed
    if dummy ~= tseries.ntimes
        % perhaps transpose?
        if settings.numberregressors == tseries.ntimes
            tseries.regressors=tseries.regressors';
            [settings.numberregressors dummy]=size(tseries.regressors);
        end
    end
else
    settings.numberregressors=0;
end
    

tseries.Y=tseries.y;



% set observations at missing epochs to NaN
if settings.multivariate
    for p=1:tseries.ntseries
        tseries.Y(p,tseries.missingepochs{p})=NaN;
        tseries.indexavailable{p}=find(~isnan(tseries.Y(p,:)));
    end
    
    
else
    if ~isempty(tseries.missingepochs)
    tseries.Y(tseries.missingepochs)=NaN;
    tseries.indexavailable=find(~isnan(tseries.Y));
    end
end

tseries.miny=min(tseries.y);
tseries.maxy=max(tseries.y);

% Check whether data is equally spaced 
% find(mod(t,1)) returns index of non-integer t values
if settings.checktime
    IndexIntegerTime=find(abs(mod(tNorm,1))<1e-3);
    if length(IndexIntegerTime)==length(tseries.t)
        settings.continuoustime=false
    else
        settings.continuoustime=true
    end
end

% See if there are data gaps in equally spaced data
if ~settings.continuoustime
   if  ~isempty(find(tseries.missingepochs))
       tseries.datagaps=true;
   end
end


% set np (number of single observations)

if ~settings.multivariate
    % univariate
    tseries.nmissingepochs=length(tseries.missingepochs);
    tseries.np=tseries.ntimes-tseries.nmissingepochs;
    % define number of epochs with actual observations
    tseries.ntimesobs=tseries.ntimes-tseries.nmissingepochs;
else
    % multivariate
    % check if there are epochs where all observations are missing
    EmptyEpochs=find(sum(isnan(tseries.y))==tseries.ntseries);
    if ~isempty(EmptyEpochs)
        % remove epochs with no observations at all
        disp(strcat('number of fully empty epochs:',num2str(length(EmptyEpochs))))
       
        tseries.ntimesobs=tseries.ntimes-length(EmptyEpochs);

        
    else
        tseries.ntimesobs=tseries.ntimes;
        
        
    end
    % total number of observations (count numbers)
    tseries.np=length(find(isnan(tseries.y)==0));
    
    if isempty(find(isnan(tseries.Y)))
        tseries.missingobs=0;
        % no missing data
    else
        tseries.missingobs=1;
        
    end
    tseries.epochsmissingobs=sum(isnan(tseries.Y))>0;
end

if settings.normalise
    % set normalized cycle periods
    settings.numbercycles=length(settings.periodscycle);

    for i=1:settings.numbercycles
       settings.lambda(i)=2*pi/settings.periodscycle(i)*tseries.period;
    end
else
    settings.numbercycles=length(settings.periodscycle);

    for i=1:settings.numbercycles
       settings.lambda(i)=2*pi/settings.periodscycle(i);
    end
end



end


