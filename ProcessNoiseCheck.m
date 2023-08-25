function [distvar,settings]=ProcessNoiseCheck(tseries,settings,distvar)
% ProcessNoiseCheck checks the completeness of the contents of distvar
%
% HOW:     [distvar]=ProcessNoiseCheck(tseries,settings,distvar)
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
%               settings.ARphi                     phi coefficient for AR
%               process (in case of settings.AR = 1)
%
%        distvar, structure with:
%           distvar.level                       variance of level term (usually zero to create smooth trend)
%           distvar.slope                       variance of slope term (zeta in Durbin&Koopman)
%           distvar.acc                         variance of acceleration
%           distvar.cycle [settings.numbercycles x 1]    variance of cycle terms
%           distvar.cyclestar [settings.numbercycles x 1]    variance of auxiliary cycle terms
%           distvar.irr                         variance of disturbance
%           distvar.AR                         variance of AR process
%
%           distvar.covirr                    covariance vector (optional)
%
% Output:
%        distvar, structure with:
%           distvar.level                       variance of level term (usually zero to create smooth trend)
%           distvar.slope                       variance of slope term (zeta in Durbin&Koopman)
%           distvar.acc                         variance of acceleration
%           distvar.cycle [settings.numbercycles x 1]    variance of cycle terms
%           distvar.cyclestar [settings.numbercycles x 1]    variance of auxiliary cycle terms
%           distvar.irr                         variance of disturbance
%           distvar.AR                         variance of AR process
%
%           distvar.covirr                    covariance vector
%
%        settings.ARphi                     phi coefficient for AR process
%
%
% Note:
%
% Taco Broerse, Utrecht University, 2018
% d.b.t.broerse@uu.nl
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
% Version 1.0 June 2018 DBT Broerse
%
%
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
%----------------------------------------------------------------------------
narginchk(3,3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do checks on process noise variances

%%%%%%%%%%%%%%%
% level
%%%%%%%%%%%%%%%

if ~isfield(distvar,'level')
    error('no level process variance set in distvar')
end

% extra check on whether single value is given
if settings.multivariate
    % Variances
    
    LengthLevel = numel(distvar.level);
    if LengthLevel == 1
        % only a single value is given, copy this for all time series
        DistVarLevel=distvar.level;
        for p=1:tseries.ntseries
            distvar.level(p)=DistVarLevel;
        end
    elseif LengthLevel == tseries.ntseries
        % ok, nothing to report
    else
        LengthLevel
        error('dimension of level process noise variance incorrect')
    end
    
    % Covariances
    
    % check whether covariances are supplied
    if isfield(distvar,'covlevel')
        [RowLevelCov,ColLevelCov] = size(distvar.covlevel);
        if RowLevelCov*ColLevelCov == 1
            % only a single value is given, copy this for all time series
            DistCovLevel=distvar.covlevel;
            distvar.covlevel=zeros(tseries.ntseries,tseries.ntseries);
            for p=1:tseries.ntseries
                for pp=1:tseries.ntseries
                    if p~=pp
                        distvar.level(p,pp)=DistCovLevel;
                    end
                end
            end
        elseif RowLevelCov == tseries.ntseries && ColLevelCov ==  tseries.ntseries
            % ok, a full matrix has been supplied
        else
            RowLevelCov
            ColLevelCov
            error('dimension of level process noise covariance incorrect')
        end
    else
        % set to zero
        distvar.covlevel=zeros(tseries.ntseries,tseries.ntseries);
    end
end

%%%%%%%%%%%%%%%
% slope
%%%%%%%%%%%%%%%
if settings.slope
    if ~isfield(distvar,'slope')
        error('no slope process variance set in distvar')
    end
    
    % Variances
    
    % extra check on whether single value is given
    if settings.multivariate
        LengthSlope = numel(distvar.slope);
        if LengthSlope == 1
            % only a single value is given, copy this for all time series
            DistVarSlope=distvar.slope;
            for p=1:tseries.ntseries
                distvar.slope(p)=DistVarSlope;
            end
        elseif LengthSlope == tseries.ntseries
            % ok
        else
            LengthSlope
            error('dimension of slope process noise variance incorrect')
        end
        
        % Covariances
        
        % check whether covariances are supplied
        if isfield(distvar,'covslope')
            [RowSlopeCov,ColSlopeCov] = size(distvar.covslope);
            if RowSlopeCov*ColSlopeCov == 1
                % only a single value is given, copy this for all time series
                DistCovSlope=distvar.covslope;
                distvar.covslope=zeros(tseries.ntseries,tseries.ntseries);
                for p=1:tseries.ntseries
                    for pp=1:tseries.ntseries
                        
                        if p~=pp
                            distvar.covslope(p,pp)=DistCovSlope;
                        end
                    end
                end
            elseif RowSlopeCov == tseries.ntseries && ColSlopeCov ==  tseries.ntseries
                % ok, a full matrix has been supplied
            else
                RowSlopeCov
                ColSlopeCov
                error('dimension of slope process noise covariance incorrect')
            end
        else
            % set to zero
            distvar.covslope=zeros(tseries.ntseries,tseries.ntseries);
        end
    end
else
    distvar.covslope=[];
    distvar.slope=[];
end

%%%%%%%%%%%%%%%
% acceleration
%%%%%%%%%%%%%%%
if settings.acc
    if ~isfield(distvar,'acc')
        error('no acceleration process variance set in distvar')
    end
    
    % Variances
    
    % extra check on whether single value is given
    if settings.multivariate
        LengthAcc = numel(distvar.acc);
        if LengthAcc == 1
            % only a single value is given, copy this for all time series
            DistVarAcc=distvar.acc;
            for p=1:tseries.ntseries
                distvar.acc(p)=DistVarAcc;
            end
        elseif LengthAcc == tseries.ntseries
            % ok
        else
            LengthAcc
            error('dimension of acceleration process noise variance incorrect')
        end
        
        % Covariances
        
        % check whether covariances are supplied
        if isfield(distvar,'covacc')
            [RowAccCov,ColAccCov] = size(distvar.covacc);
            if RowAccCov*ColAccCov == 1
                % only a single value is given, copy this for all time series
                DistCovAcc=distvar.covacc;
                distvar.covacc=zeros(tseries.ntseries,tseries.ntseries);
                for p=1:tseries.ntseries
                    for pp=1:tseries.ntseries
                        if p~=pp
                            distvar.covacc(p,pp)=DistCovAcc;
                        end
                    end
                end
            elseif RowAccCov == tseries.ntseries && ColAccCov ==  tseries.ntseries
                % ok, a full matrix has been supplied
            else
                RowAccCov
                ColAccCov
                error('dimension of acceleration process noise covariance incorrect')
            end
        else
            % set to zero
            distvar.covacc=zeros(tseries.ntseries,tseries.ntseries);
        end
    end
else
    distvar.covacc=[];
    distvar.acc=[];
end

%%%%%%%%%%%%%%%
% cycles
%%%%%%%%%%%%%%%

% cycles are currently treated univariately, with no covariances between
% the cycles of the different timeseries

if settings.cycle
    if ~isfield(distvar,'cycle')
        error('no cycle process variance set in distvar')
    end
    
    
    % Variances
    
    % extra check on whether single value is given
    if settings.multivariate
        % loop on cycles
        cyclecopy=distvar.cycle;
        if isfield(distvar,'covcyclestar')
            covcyclecopy=distvar.covcyclestar;
        else
            covcyclecopy=[];
        end
        
        for ii=1:settings.numbercycles
            LengthCycle = numel(cyclecopy(:,ii));
            if LengthCycle == 1
                % only a single value is given, copy this for all time series
                DistVarCycle=cyclecopy(1,ii);
                
                for p=1:tseries.ntseries
                    distvar.cycle(p,ii)=DistVarCycle;
                end
            elseif LengthCycle == tseries.ntseries
                % ok
            else
                LengthCycle
                error('dimension of cycle process noise variance incorrect')
            end
            
            % check whether covariances are supplied between
            % cycle and cyclestar
            if ~isfield(distvar,'covcyclestar') || isempty(covcyclecopy)
                for p=1:tseries.ntseries
                    distvar.covcyclestar(p,ii)=0;
                end
            else
                [RowCycleCovStar,ColCycleCovStar] = size(covcyclecopy);
                % row should coincide with number of timeseries or be one
                % column with number of cycles
                if ColCycleCovStar ~= settings.numbercycles
                    error('number of columns in distvar.covcyclestar has to coincide with the number of cycle terms')
                end
                
                % row should coincide with number of timeseries or be one
                % column with number of cycles
                
                if RowCycleCovStar == 1
                    % only a single value is given, copy this for all time series
                    DistVarCycleCovStar=covcyclecopy(1,ii);
                    
                    for p=1:tseries.ntseries
                        
                        distvar.covcyclestar(p,ii)=DistVarCycleCovStar;
                    end
                elseif RowCycleCovStar == tseries.ntseries
                    % ok
                else
                    RowCycleCovStar
                    error('dimension of covariance between cycle and cycle star process noise variance incorrect')
                end
            end
        end
    end
    
    % check if a separate cycle star variance has been supplied
    if ~isfield(distvar,'cyclestar')
        for p=1:tseries.ntseries
            for ii=1:settings.numbercycles
                distvar.cyclestar(p,ii)=distvar.cycle(p,ii);
            end
        end
    else
        if settings.multivariate
            % loop on cycle star
            cyclestarcopy=distvar.cyclestar;
            
            
            for ii=1:settings.numbercycles
                LengthCycle = numel(cyclestarcopy(:,ii));
                if LengthCycle == 1
                    % only a single value is given, copy this for all time series
                    DistVarCycle=cyclestarcopy(1,ii);
                    
                    for p=1:tseries.ntseries
                        distvar.cyclestar(p,ii)=DistVarCycle;
                    end
                elseif LengthCycle == tseries.ntseries
                    % ok
                else
                    LengthCycle
                    error('dimension of cycle star process noise variance incorrect')
                end
            end
        end
    end
    % check if there are already covariances set for cycle and cyclestar
    if ~isfield(distvar,'covcyclestar') || isempty(distvar.covcyclestar)
        for ii=1:settings.numbercycles
            distvar.covcyclestar(1,ii)=0;
        end
    end
else
    distvar.covcyclestar=[];
    distvar.cycle=[];
end



%%%%%%%%%%%%%%%
% AR
%%%%%%%%%%%%%%%

if settings.AR
    if ~isfield(distvar,'AR')
        error('no AR process variance set in distvar')
    end
    
    % Variances
    
    % extra check on whether single value is given
    if settings.multivariate
        LengthAR = numel(distvar.AR);
        if LengthAR == 1
            % only a single value is given, copy this for all time series
            DistVarAR=distvar.AR;
            for p=1:tseries.ntseries
                distvar.AR(p)=DistVarAR;
            end
        elseif LengthAR == tseries.ntseries
            % ok
        else
            LengthAR
            error('dimension of acceleration process noise variance incorrect')
        end
        
        % Covariances
        
        % check whether covariances are supplied
        if isfield(distvar,'covAR')
            [RowARCov,ColARCov] = size(distvar.covAR);
            if RowARCov*ColARCov == 1
                % only a single value is given, copy this for all time series
                DistCovAR=distvar.covAR;
                distvar.covAR=zeros(tseries.ntseries,tseries.ntseries);
                for p=1:tseries.ntseries
                    for pp=1:tseries.ntseries
                        if p~=pp
                            distvar.covAR(p,pp)=DistCovAR;
                        end
                    end
                end
            elseif RowARCov == tseries.ntseries && ColARCov ==  tseries.ntseries
                % ok, a full matrix has been supplied
            else
                RowARCov
                ColARCov
                error('dimension of AR process noise covariance incorrect')
            end
        else
            % set to zero
            distvar.covAR=zeros(tseries.ntseries,tseries.ntseries);
        end
        
        % phi values
        LengthPhi = length(settings.ARphi);
        if LengthPhi == 1
            ARphi=settings.ARphi;
            % make settings.ARphi as long as number of time series
            settings.ARphi([1:tseries.ntseries],1)=ARphi;
        elseif LengthPhi==tseries.ntseries
            % ok
            [a,b]=size(settings.ARphi);
            if b==settings.ARorder
                % ok
            elseif a==settings.ARorder
                % switch dimensions
                settings.ARphi=settings.ARphi';
            else
                
                error('number of AR phi parameters not equal to AR order')
            end
        else
            LengthPhi
            error('length of settings.ARphi incorrect, this should either be 1 or equal to the number of time series')
        end
    end
else
    distvar.AR=[];
    distvar.covAR=[];
end

%%%%%%%%%%%%%%%
% irregular term
%%%%%%%%%%%%%%%

if ~isfield(distvar,'irr')
    error('no irregular process variance set in distvar')
end

% Variances

% extra check on whether single value is given
if settings.multivariate
    LengthIrr = numel(distvar.irr);
    if LengthIrr == 1
        % only a single value is given, copy this for all time series
        DistVarIrr=distvar.irr;
        for p=1:tseries.ntseries
            distvar.irr(p)=DistVarIrr;
        end
    elseif LengthIrr == tseries.ntseries
        % ok
    else
        LengthIrr
        error('dimension of irregular process noise variance incorrect')
    end
    
    % Covariances
    
    % check whether covariances are supplied
    if isfield(distvar,'covirr')
        [RowAccIrr,ColAccIrr] = size(distvar.covirr);
        if RowAccIrr*ColAccIrr == 1
            % only a single value is given, copy this for all time series
            DistCovIrr=distvar.covirr;
            distvar.covirr=zeros(tseries.ntseries,tseries.ntseries);
            for p=1:tseries.ntseries
                for pp=p+1:tseries.ntseries
                    %                     p
                    %                     pp
                    %                     if p~=pp
                    distvar.covirr(pp,p)=DistCovIrr;
                    distvar.covirr(p,pp)=DistCovIrr;
                    %                     end
                end
            end
            
        elseif RowAccIrr == tseries.ntseries && ColAccIrr ==  tseries.ntseries
            % ok, a full matrix has been supplied
        else
            RowAccIrr
            ColAccIrr
            error('dimension of irregular process noise covariance incorrect')
        end
    else
        % set to zero
        distvar.covirr=zeros(tseries.ntseries,tseries.ntseries);
    end
else
    % set to zero
    distvar.covirr=zeros(tseries.ntseries,tseries.ntseries);
end

% variance and covariance to a single matrix
distvar.varcovirr=diag(distvar.irr)+distvar.covirr;


end


