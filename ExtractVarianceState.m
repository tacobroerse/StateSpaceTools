function [var]=ExtractVariance(VarMatrix,tseries,settings,model,denormalise)
% ExtractVarianceState.m extracts the corresponding variances for each state
% space component from a combined process variance matrix.
%
%
% HOW:      [var]=ExtractVariance(VarMatrix,tseries,settings,model)
%           [var]=ExtractVariance(VarMatrix,tseries,settings,model,denormalise)
%
% Input:    VarMatrix                       process variance matrix
%
%      tseries, structure of which this function uses:
%           tseries.period                  time series sampling periodicity
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
%           model                           contains indices of components
%
%
%
% Output:   var.level [p x n]                variance of level term (usually zero to create smooth trend)
%           var.slope [p x n]                variance of slope term (zeta in Durbin&Koopman)
%           var.acc   [p x n]                variance of acceleration
%           var.cycle [settings.numbercycles x 2 x n] variance of cycle terms
%           var.delta [n_interventions]  variance of step interventions
%           var.regressors  [n_regressors]     variance of regression coefficients
%           var.AR [p x n]                   variance of AR
%            where p is the number of time series
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
% Version 1.0 june 2014 DBT Broerse
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
% Version 2.1 June 2018 DBT Broerse
% - extended to multivariate models
% - updated use of model.indices
%
%----------------------------------------------------------------------------
% remarks: Time dependence in H and Q is not properly accounted for
%----------------------------------------------------------------------------
%disp('introduce time variance in extract variance')
narginchk(4,5);

% set defaults
if isempty(settings.type);settings.type='Normal_Slope_and_Acc';end
if nargin==4;denormalise=false;end

%[Dimension,dummy,ntimes]=size(VarMatrix);

ntimes=tseries.ntimes;

% set up matrices
for p=1:tseries.ntseries
    var.level{p}=zeros(ntimes);var.slope=[];var.cycle=[];var.acc=[];var.delta=[];var.regressors=[];
    if settings.slope
        var.slope{p}=zeros(1,ntimes);
    end
    if settings.acc
        var.acc{p}=zeros(1,ntimes);
    end
    if settings.cycle
        var.cycle{p}=zeros(settings.numbercycles,ntimes);
    end
    if settings.AR
        var.AR{p}=zeros(1,ntimes);
    end
end


%     if strcmp(settings.type,'Lag_Operator_Both_Slope_and_Acc')
%         %% Determine dimension of matrices
%         % trend adds 2 dimensions
%         % acceleration 3
%         % each cycle 2
%
%         % trend
%         var.slope(1,:)=VarMatrix(1,1,:);
%         NextDim=3;
%         % quadratic term
%         if settings.acc
%            var.acc(1,:)=VarMatrix(3,3,:);
%            NextDim=NextDim+3;
%         end
%
%         % cycle terms
%         if settings.cycle
%             CycleDim=NextDim;
%             for i=1:settings.numbercycles
%                 var.cycle(1,i,:)=VarMatrix(CycleDim+(i-1)*2,CycleDim+(i-1)*2,:);
%                 var.cycle(2,i,:)=VarMatrix(CycleDim+(i-1)*2+1,CycleDim+(i-1)*2+1,:);
%                 %Q(CycleDim+(i-1)*2+1,CycleDim+(i-1)*2+1)=var.cycle(i);
%             end
%         end
%
%
%     if strcmp(settings.type,'Lag_Operator_either_Slope_or_Acc')
%
%         if settings.acc
%             % quadratic term
%            var.acc(1,:)=VarMatrix(1,1,:);
%            NextDim=4;
%         else
%             % trend
%            var.slope(1,:)=VarMatrix(1,1,:);
%            NextDim=3;
%         end
%
%         % cycle terms
%         if settings.cycle
%             for i=1:settings.numbercycles
%                 var.cycle(1,i,:)=VarMatrix(CycleDim+(i-1)*2,CycleDim+(i-1)*2,:);
%                 var.cycle(2,i,:)=VarMatrix(CycleDim+(i-1)*2+1,CycleDim+(i-1)*2+1,:);
%             end
%         end
%
%
%
if strcmp(settings.type,'Normal_Slope_and_Acc')
    
    
    for p=1:tseries.ntseries
        var.level{p}=squeeze(VarMatrix(model.indexlevel{p},model.indexlevel{p},:))';
      
        
        if settings.slope
            var.slope{p}=squeeze(VarMatrix(model.indexslope{p},model.indexslope{p},:))';
            if settings.acc
                var.acc{p}=squeeze(VarMatrix(model.indexacc{p},model.indexacc{p},:))';
            end
        end
        
        % cycle terms
        if settings.cycle
            for i=1:settings.numbercycles
                var.cycle{p}(i,:)=squeeze(VarMatrix(model.indexcycles{p}(i*2-1),model.indexcycles{p}(i*2-1),:))';
                %    var.cyclestar{p}(i,:)=VarMatrix(CycleDim+(i-1)*2+1,CycleDim+(i-1)*2+1,:);
            end
        end
        
        if settings.regression
            for i=1:settings.numberregressors
                var.regressors{p}(i)=squeeze(VarMatrix(model.indexregressors{p}(i),model.indexregressors{p}(i),1))';
            end
        end
        
        if settings.intervention
            for i=1:settings.numberinterventions
                
                var.delta{p}(i)=squeeze(VarMatrix(model.indexinterventions{p}(i),model.indexinterventions{p}(i),1))';
            end
        end
        
        if settings.AR
            var.AR{p}=squeeze(VarMatrix(model.indexAR{p},model.indexAR{p},:))';
        end
    end
    
else
    disp('invalid choice for state space formulation')
    settings.type
    return
end

if settings.normalise
    % divide all variances by time period for denormalisation
    if denormalise
        
        disp('denormalising of state error variances in ExtractVarianceState')
        for p=1:tseries.ntseries
            if settings.slope
                var.slope{p}=var.slope{p}/tseries.period^2;
                if settings.acc
                    % var.acc=var.acc/tseries.period^2;
                    var.acc{p}=var.acc{p}/tseries.period^4;
                    disp('check acceleration variance denormalisation')
                end
                
                % other variances have no time dependency
                
                %                 if settings.cycle
                %
                %                     for i=1:settings.numbercycles
                %                         %  var.cycle(1,i,:)=var.cycle(1,i,:)/tseries.period^2;
                %                         % var.cycle(2,i,:)=var.cycle(2,i,:)/tseries.period^2;
                %                     end
                %                 end
                %
                %                 if settings.regression
                %                     for i=1:settings.numberregressors
                %
                %                         %  var.regressors(i)=var.regressors(i)/tseries.period^2;
                %
                %                     end
                %                 end
                %
                %                 if settings.intervention
                %                     for i=1:settings.numberinterventions
                %
                %                         % var.delta(i)=var.delta(i)/tseries.period^2;
                %
                %                     end
                %                 end
            end
        end
    end
end





end