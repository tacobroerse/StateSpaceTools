function [var]=ExtractVarianceDisturbance(VarEta,VarEpsilon,model,tseries,settings)
% ExtractVariance.m extracts the corresponding variances for each state
% space component from a combined process variance matrix.
%
%
% HOW:      [var]=ExtractVarianceDisturbance(VarEta,VarEtaIrregular,model,tseries,settings)
%
% Input:    VarEta                       process covariance matrix
%
%           VarEpsilon                   irregular covariance matrix
%
%      model                                contain index description of components
%
%      tseries, structure of which this function uses:
%           tseries.period                  time series sampling periodicity
%           tseries.dt                      time step
%           ntimes                  number of epochs
%
%      settings, structure of which this function uses:
%           settings.slope                  determines whether a slope will be estimated
%           settings.acc                    determines whether accelerations will be estimated
%           settings.cycle                  determines whether cycle components will be estimated
%           settings.regression             determines whether regression is used (using user supplied regressors)
%           settings.intervention           determines whether step intervention effects will be estimated
%           settings.numbercycles           number of cycle terms
%           settings.numberregressors       number of regressors
%           settings.numberinterventions    number of interventions
%           settings.AR                     determines whether AR components will be estimated
%
%           settings.type                   type of formulation of the state space equations
%                                           - 'Lag_Operator_Both_Slope_and_Trend'
%                                           using a lag operator and separate components
%                                           for slope and acceleration
%                                           - 'Lag_Operator_either_Slope_or_Trend'
%                                           using a lag operator and combined components
%                                           for slope and acceleration
%                                           - 'Normal_Slope_and_Acc' default
%
%           denormalise                     denormalise trend and acceleration variance
%
%
%
% Output:   var.level [p]                variance of level term (usually zero to create smooth trend)
%           var.slope [p]                variance of slope term (zeta in Durbin&Koopman)
%           var.acc   [p]                variance of acceleration
%           var.cycle [p x settings.numbercycles] variance of cycle terms
%           var.cyclestar [p x settings.numbercycles] variance of cycle terms
%           var.covcyclestar [p x settings.numbercycles] covariance of cycle and cycle star terms
%           coefficients
%           var.AR [p]                     variance of AR process
%           var.irr [p]                  variance of irregular term
%           var.eta [n]                  variance of state vector disturbances
%           var.covirr [p x p]           covariance of irregular term (for multivariate models)
%           var.varcorirr [p x p]        variance-covariance of irregular term (for multivariate models)
%           var.covlevel [p x p]         covariance of level term (for multivariate models)
%           var.covslope [p x p]         covariance of slope term (for multivariate models)
%           var.covacc   [p x p]         covariance of acceleration (for multivariate models)
%           var.covAR  [p x p]           covariance of AR process (for multivariate models)
%           var.varcoreta [n x n]        variance-covariance of disturbances (for multivariate models)
%
%
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
% added: allows now for extra dimension (time)
%
% Version 1.2 August 2014 DBT Broerse
% added: option for denormalizing derivatives of the trend
%
% Version 1.3 September 2014 DBT Broerse
% - added regression
% - added pulse intervention
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
% - continuous time and integer time functions merged
%
% Version 3.0 April 2014 DBT Broerse
% - introduced time component in matrices
% - phased out state space definitions based on lag operator
% - separate subroutine for disturbance noise variance and state variance
%
% Version 3.1 June 2018 DBT Broerse and T Frederikse
% - added covariance for cycle & cycle star
% - added variance for AR
%
% Version 3.2 June 2018 DBT Broerse
% - extension to multivariate models
% - addition of extraction of variance of irregular term (epsilon)
%
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------
narginchk(5,6);


% set defaults
if isempty(settings.type);settings.type='Normal_Slope_and_Acc';end
if nargin==5;denormalise=false;end

%[Dimension,dummy,ntimes]=size(VarEta);

dt=tseries.dt;

% set up matrices
% variances
var.level=[];var.slope=[];var.cycle=[];var.cyclestar=[];var.acc=[];var.AR=[];var.irr=[];
% variance vector
var.eta=[];
% covariances
var.covlevel=[];var.covslope=[];var.covcyclestar=[];var.covacc=[];var.covAR=[];var.covirr=[];
% correlations
var.corrlevel=[];var.corrslope=[];var.corrcyclestar=[];var.corracc=[];var.corrAR=[];var.corrirr=[];
% full variance covariance matrices
var.varcovirr=[];var.varcoveta=[];
if strcmp(settings.type,'Normal_Slope_and_Acc')
    for p=1:tseries.ntseries
        
        
        
        % irregular term
        var.irr(p)=VarEpsilon(p,p);
        
        % state variables:
        
        % variance
        var.level(p)=VarEta(model.indexlevel{p},model.indexlevel{p});
        if settings.slope
            var.slope(p)=VarEta(model.indexslope{p},model.indexslope{p});
            if settings.acc
                var.acc(p)=VarEta(model.indexacc{p},model.indexacc{p});
            end
        end
        
        % cycle terms
        if settings.cycle
            for i=1:settings.numbercycles
                var.cycle(p,i)=VarEta(model.indexcycles{p}(2*i-1),model.indexcycles{p}(2*i-1));
                var.cyclestar(p,i)=VarEta(model.indexcycles{p}(2*i),model.indexcycles{p}(2*i));
                % covariance between cycle and cyclestar
                if settings.distvarcovariancecycles
                    var.covcyclestar(p,i)=VarEta(model.indexcycles{p}(2*i-1),model.indexcycles{p}(2*i));
                    % correlation
                    var.corrcyclestar(p,i)=var.covcyclestar(p,i)/sqrt(var.cycle(p,i)*var.cyclestar(p,i));
                else
                    var.covcyclestar(p,i)=0;
                    var.corrcyclestar(p,i)=0;
                end
                
            end
        end
        
        
        
        if settings.AR
            %  for i=1:settings.ARorder
            var.AR(p)=VarEta(model.indexAR{p},model.indexAR{p});
            %  end
        end
        
        if settings.multivariate
            % add covariances as well
            for pp=p+1:tseries.ntseries
                
                if settings.multivariatedistcovarianceirr
                    % irregular term
                    var.covirr(p,pp)=VarEpsilon(p,pp);
                    var.covirr(pp,p)=VarEpsilon(pp,p);
                    
                    % correlation
                    var.corrirr(p,pp)=VarEpsilon(p,pp)/sqrt(VarEpsilon(p,p)*VarEpsilon(pp,pp));
                    var.corrirr(pp,p)=var.corrirr(p,pp);
                else
                    var.covirr(p,pp)=0;
                    var.covirr(pp,p)=0;
                    var.corrirr(p,pp)=0;
                    var.corrirr(pp,p)=0;
                end
                
                if settings.multivariatedistcovariancelevel
                    % level covariances
                    var.covlevel(p,pp)=VarEta(model.indexlevel{p},model.indexlevel{pp});
                    var.covlevel(pp,p)=VarEta(model.indexlevel{pp},model.indexlevel{p});
                    % correlation
                    
                    var.corrlevel(p,pp)=VarEta(p,pp)/sqrt(VarEta(model.indexlevel{p},model.indexlevel{p})*...
                        VarEta(model.indexlevel{pp},model.indexlevel{pp}));
                    var.corrlevel(pp,p)=var.corrlevel(p,pp);
                else
                    % level covariances
                    var.covlevel(p,pp)=0;
                    var.covlevel(pp,p)=0;
                    % correlation
                    
                    var.corrlevel(p,pp)=0;
                    var.corrlevel(pp,p)=0;
                end
                % slope
                if settings.slope
                    if settings.multivariatedistcovarianceslope
                        var.covslope(p,pp)=VarEta(model.indexslope{p},model.indexslope{pp});
                        var.covslope(pp,p)=VarEta(model.indexslope{pp},model.indexslope{p});
                        % correlation
                        
                        var.corrslope(p,pp)=var.covslope(p,pp)/sqrt(VarEta(model.indexslope{p},model.indexslope{p})*...
                            VarEta(model.indexslope{pp},model.indexslope{pp}));
                        var.corrslope(pp,p)=var.corrslope(p,pp);
                    else
                        var.covslope(p,pp)=0;
                        var.covslope(pp,p)=0;
                        % correlation
                        
                        var.corrslope(p,pp)=0;
                        var.corrslope(pp,p)=0;
                    end
                end
                % acceleration
                if settings.acc
                    if settings.multivariatedistcovarianceacc
                        var.covacc(p,pp)=VarEta(model.indexacc{p},model.indexacc{pp});
                        var.covacc(pp,p)=VarEta(model.indexacc{pp},model.indexacc{p});
                        % correlation
                        
                        var.corracc(p,pp)=var.covacc(p,pp)/sqrt(VarEta(model.indexacc{p},model.indexacc{p})*...
                            VarEta(model.indexacc{pp},model.indexacc{pp}));
                        var.corracc(pp,p)=var.corracc(p,pp);
                    else
                        var.covacc(p,pp)=0;
                        var.covacc(pp,p)=0;
                        % correlation
                        
                        var.corracc(p,pp)=0;
                        var.corracc(pp,p)=0;
                    end
                    
                end
                % cycle terms are currently not included for multivariate
                % correlation
                
                
                
                if settings.AR
                    if settings.multivariatedistcovarianceAR
                        % for i=1:settings.ARorder
                        var.covAR(p,pp)=VarEta(model.indexAR{p},model.indexAR{pp});
                        var.covAR(pp,p)=VarEta(model.indexAR{pp},model.indexAR{p});
                        %  end
                        % correlation
                        
                        var.corrAR(p,pp)=var.covAR(p,pp)/sqrt(VarEta(model.indexAR{p},model.indexAR{p})*...
                            VarEta(model.indexAR{pp},model.indexAR{pp}));
                        var.corrAR(pp,p)=var.corrAR(p,pp);
                    else
                        % for i=1:settings.ARorder
                        var.covAR(p,pp)=0;
                        var.covAR(pp,p)=0;
                        %  end
                        % correlation
                        
                        var.corrAR(p,pp)=0;
                        var.corrAR(pp,p)=0;
                    end
                end
            end
            
            
            
            
        end
    end
    % add variances and covariances
    if settings.multivariate
        if ~isempty(var.covirr)
            var.varcovirr=diag(var.irr)+var.covirr;
        else
            var.varcovirr=diag(var.irr);
        end
    else
        var.varcovirr=(var.irr);
    end
    
    % the full variance-covariance matrix is simply VarEta
    var.varcoveta=VarEta;
    % variance is the diagonal
    var.eta=diag(VarEta);
    
else
    disp('invalid choice for state space formulation')
    settings.type
    return
end

%     if settings.normalise
%         disp('hop je sokkel mop')
%
%         hop
%         % divide all variances by time period for denormalisation
%         if denormalise
%             disp('denormalising')
%             if settings.slope
%                 var.slope=var.slope/tseries.period^2;
%                 if settings.acc
%                     var.acc=var.acc/tseries.period^2;
%                 end
%                 if settings.cycle
%
%                     for i=1:settings.numbercycles
%                         var.cycle(p,i)=var.cycle(p,i)/tseries.period^2;
%                         var.cyclestar(p,i)=var.cyclestar(p,i)/tseries.period^2;
%                         var.covcyclestar(p,i)=var.covcyclestar(p,i)/tseries.period^2;
%                     end
%                 end
%
% %                 if settings.regression
% %                     for i=1:settings.numberregressors
% %
% %                         var.regressors(p,i)=var.regressors(i)/tseries.period^2;
% %
% %                     end
% %                 end
% %
% %                 if settings.intervention
% %                     for i=1:settings.numberinterventions
% %
% %                         var.delta(i)=var.delta(i)/tseries.period^2;
% %
% %                     end
% %                 end
%
%                 if settings.AR
%                     for i=1:settings.ARorder
%
%                         var.AR(p,i)=var.AR(p,i)/tseries.period^2;
%
%                     end
%                 end
%             end
%         end
%     end




end