function [tseries]=FindMissingEpochs(tseries)
% FindMissingEpochs.m finds missing epochs for equally spaced data.
% HOW:     [time,y,MissingEpochs]=FindMissingEpochs(time,y)
%   
% Input:    tseries structure which should contain:
%               tseries.time [epochs x 1]           time vector
%               tseries.y [epochs x 1]              observation vector           
%
% Output:   
%  
%
% Note: 
%
% Taco Broerse, Delft University of Technology, 2014
% d.b.time.broerse@tudelft.nl
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
% Version 2.0 October 2014 DBT Broerse
% - changed I/O
%
%----------------------------------------------------------------------------
% remarks: 
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
%----------------------------------------------------------------------------
narginchk(1,1)

% extract from structure tseries
y=tseries.y;
time=tseries.time;

% defaults
NMissingEpochs=0;
eps=1e-3;

% number of epochs
ntimes=length(time);
% check time steps
dt=diff(time);

% take minimum time step as default
Defaultdt = min(dt);


% find time steps larger than default dt

Gaps=find(abs(dt-Defaultdt)>1e-3*Defaultdt);
if ~isempty(Gaps)
    GapsLeft=true;
    % where to start looking for gaps
    StartSearch=2;

    while(GapsLeft)

        % start searching for gaps
        for i=StartSearch:ntimes-1
%             
%             if (NMissingEpochs>20)
%                 disp('we stoppen even')
%                 return
%             end
            if (abs(dt(i)-Defaultdt)>eps*Defaultdt)
                
                StartSearch=i-1;
                NMissingEpochs=NMissingEpochs+1;
                % add a time step
                time(ntimes+1)=time(i)+Defaultdt;
                y(ntimes+1)=NaN;
                % sort
                [time,index]=sort(time,1);
                y=y(index);
                % determine time steps again
                dt=diff(time);
                % new number of epochs
                ntimes=length(time);
                % quit loop and start again
              
                break
                
            end
        end

        % check again for gaps
        Gaps=find(abs(dt-Defaultdt)>1e-3*Defaultdt);
        if isempty(Gaps)
            GapsLeft=false;
        end
        
    end

    tseries.missingepochs=find(isnan(y));
else
    tseries.missingepochs=[];

end

% save in structure

tseries.time=time;
tseries.y=y;
tseries.ntimes=ntimes;


end

