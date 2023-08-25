%function [fit,trend,trend_no_acc,cycles,phases,residuals,residualRMS,bias,slope,acceleration] = LeastSquaresTrend(time,y,Op_Acc,Op_Cycle,NumberCycles,lambda)
function [estlsq]=LeastSquaresTrend(tseries,settings,p)

% LeastSquaresTrend calculates deterministic fit 
% 
%
% HOW:     [estlsq]=LeastSquaresTrend(tseries,settings)
%
% Input:    structure tseries, of which is used:
%               tseries.time [epochs x 1]           time vector
%               tseries.y                            observation vector
%
%           structure settings, which should include:
%               settings.acclsq                   determines whether accelerations will be estimated
%               settings.cycle                    determines whether cycle components will be estimated
%               settings.periodscycle             periods of the cycles
%
%
%
% Output:   structure estlsq, which includes
%               estlsq.fit                         least squares fit
%               estlsq.trend   [1 x epochs]        trend fit
%               estlsq.cycles  [NumberCycles x epochs] cycle fits
%               estlsq.phases                      phase of cycle fit
%               estlsq.residuals                   residual
%               estlsq.residualRMS                 root mean square of residuals
%               estlsq.bias
%               estlsq.slope                       
%               estlsq.acceleration
% 
%
% Note:    
%
% Taco Broerse, Delft University of Technology, 2014
% d.b.t.broerse@tudelft.nl
%
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
% - added output
% - subract reference epoch from time
%
% Version 2.0 October 2014 DBT Broerse
% - changed I/O, using structures
%
% Version 3.0 May 2015 DBT Broerse
% - added error covariance matrix
%
%----------------------------------------------------------------------------
% remarks:
%----------------------------------------------------------------------------
if nargin==2
    p=1;
end

trend=[];cycles=[];phases=[];residuals=[];bias=[];slope=[];acceleration=[];

% getting elements from structures
time=tseries.time;
if iscell(tseries.y)
    y=tseries.y{p}';
else
y=tseries.y;
end

% getting information on cycles
NumberCycles=length(settings.periodscycle);
lambda=2*pi./settings.periodscycle;

% % removing NaN from time series
% IndexAvailable=find(~isnan(y))
% time=time(IndexAvailable);
% y=y(IndexAvailable);

epochs=length(time);

%disp('no data weighting yet')

%% H: disturbance variance
%H = SDObs.^2;

%% reference epoch
t0=time(1);
time=time-t0;
%% A: normal matrix
A=[];

A(:,1) = ones(epochs,1); % to incorporate bias
A(:,2) = time';
NextDim=3;
%number of parameters
nr_param=2;
if settings.acclsq
    nr_param=nr_param+1;
    A(:,NextDim) = time.*time; %(eq.1 and eq.3) to incorporate quadratic term mu2_t+1
    NextDim=NextDim+1;
end

if settings.cycle
   CycleDim=NextDim;
   for i=1:NumberCycles
      nr_param=nr_param+2;
      NextDim=NextDim+2;
      A(:,CycleDim+(i-1)*2) = cos(lambda(i)*time'); % cosine term
      A(:,CycleDim+(i-1)*2+1) = sin(lambda(i)*time'); % sine term     
   end
end

%% least square inversion
NM = (A'*A)\A';


% normal equation

x = NM*y;

%% determine fit and split fit in several components

fit=A*x;

%trend
ATrend=zeros(size(A));
ATrend(:,1:2)=A(:,1:2);
trend=ATrend*x;  
bias=x(1);
slope=x(2);

%acceleration
if settings.acclsq
   ATrend(:,3)=A(:,3);
   trend=ATrend*x;
   acceleration=x(3);
   
end

%fit of cycle terms
if settings.cycle
   cycles=zeros(NumberCycles,epochs);
   for i=1:NumberCycles
     ACycle=zeros(size(A));
     ACycle(:,CycleDim+(i-1)*2:CycleDim+(i-1)*2+1)=A(:,CycleDim+(i-1)*2:CycleDim+(i-1)*2+1);
     cycles(i,:)=ACycle*x;
     
   end
end

residuals=y-fit;
residualRMS=sqrt(1/epochs*sum(residuals.^2));

%% variance covariance
Exx  = (A'*A)\eye(size(A'*A));

% aposteriori variance estimate: Error^2/(no of obs - no of param)
var_est = residuals'*residuals/(epochs - nr_param);
% covariance matrix for parameters
Exx_scal = var_est*Exx;


%% set up output structure

estlsq.fit = fit;                       
estlsq.trend = trend;
estlsq.cycles = cycles; 
estlsq.phases = phases;
estlsq.residuals = residuals;
estlsq.residualRMS = residualRMS;
estlsq.bias = bias;
estlsq.slope = slope;                     
estlsq.acceleration = acceleration;
estlsq.nr_param=nr_param;
estlsq.cov=Exx_scal;

end