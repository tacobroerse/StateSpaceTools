function [f,Ypsd] = SpectralAnalysis(y,time,Frequency)
% Spectral Analysis 
% author: D.B.T. Broerse, based on a matlab example

%
% HOW:     [f,Ypsd] = SpectralAnalysis(y,time)
%
% Input:    time [epochs x 1]       time vector
%           y [epochs x 1]          time series 
%           Frequency               measurement frequency
%
% Output:   f                       absolute frequencies
%           Ypsd                    power spectrum density of time series
%
% Note: this version requires equally spaced and normalized time steps, 
%       i.e. t=1,2,....N
%
% Taco Broerse, Delft University of Technology, 2014
% d.b.t.broerse@tudelft.nl
%
%----------------------------------------------------------------------------
% uses: none
% 
%----------------------------------------------------------------------------
% revision history
%
% Version 1.0 june 2014 DBT Broerse
%
%----------------------------------------------------------------------------
% remarks
%----------------------------------------------------------------------------
%
% do a check on input time
ntimes=length(time);
Period=(time(end)-time(1))/(ntimes-1);
% Normalized time (integers in case of equally spaced time)
t=(time-time(1))/Period;
% is time indeed equally spaced
if find(abs(mod(t,1))<1e-10)
    ContinuousTime=false;
else
    ContinuousTime=true;
    disp('error, time is not equally spaced')
    return
end

%% spectral analysis 

NFFT = 2^nextpow2(ntimes); % Next power of 2 from length of t
Yspec = fft(y,NFFT)/ntimes;
%Yspec = fft(y)/ntimes;
f = Frequency/2*linspace(0,1,NFFT/2+1);
%f = Frequency/2*linspace(0,1,ntimes/2+1);
Ypsd=2*abs(Yspec(1:NFFT/2+1));
%Ypsd=2*abs(Yspec(1:ntimes/2+1));


end