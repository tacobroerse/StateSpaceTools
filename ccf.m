function [ta,lags] = ccf(y1,y2,p)
% ACF - Compute Crosscorrelations Through p Lags
% >> myacf = ccf(y1,y2,p) 
%
% Inputs:
% y1 - series to compute ccf for, nx1 column vector
% y2 - series to compute ccf for, nx1 column vector
% p - total number of lags, 1x1 integer
%
% Output:
% myacf - p+1x1 vector containing crosscorrelations
% lags - p+1x1 vector containing lags, starting at zero    
%
%
% A bar graph of the crosscorrelations is also produced, with
% rejection region bands for testing individual crosscorrelations = 0.
%
%
% compared to original acf function, scaling has been changed to prevent
% damping at increasing lag 
%
% Example:
% >> ccf(randn(100,1),randn(100,1), 10)
%


% --------------------------
% USER INPUT CHECKS
% --------------------------

[n1, n2] = size(y1) ;
if n2 ~=1
    error('Input series y1 must be an nx1 column vector')
end

[m1, m2] = size(y2) ;
if m2 ~=1
    error('Input series y2 must be an nx1 column vector')
end

if m1 ~= n1
   error('Input series y1 and y2 should have same length') 
end

[a1, a2] = size(p) ;
if ~((a1==1 & a2==1) & (p<n1))
    error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
end



% -------------
% BEGIN CODE
% -------------

ta = zeros(p+1,1) ;
global N 
N = max(size(y1)) ;
global y1bar 
y1bar = mean(y1); 
global y2bar
y2bar = mean(y2);

lags=[0:p];
% Collect CCFs at each lag i
for i = 0:p
   ta(i+1) = acf_k(y1,y2,i) ; 
end

% Plot ACF
% Plot rejection region lines for test of individual autocorrelations
% H_0: rho(tau) = 0 at alpha=.05
bar(ta)
line([0 p+.5], (2.32)*(1/sqrt(N))*ones(1,2))
line([0 p+.5], (-2.32)*(1/sqrt(N))*ones(1,2))

% Some figure properties
line_hi = (2.32)*(1/sqrt(N))+.05;
line_lo = -(2.32)*(1/sqrt(N))-.05;
bar_hi = max(ta)+.1 ;
bar_lo = -max(ta)-.1 ;

if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
    axis([0 p+.60 line_lo line_hi])
else
    axis([0 p+.60 bar_lo bar_hi])
end
title({' ','Sample Autocorrelations',' '})
xlabel('Lag Length')
set(gca,'YTick',[-1:.20:1])
% set number of lag labels shown
if (p<28 & p>4)
    set(gca,'XTick',floor(linspace(1,p,4)))
elseif (p>=28)
    set(gca,'XTick',floor(linspace(1,p,8)))
end
set(gca,'TickLength',[0 0])




% ---------------
% SUB FUNCTION
% ---------------
function ta2 = acf_k(y1,y2,k)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
% 
global y1bar y2bar
global N
cross_sum = zeros(N-k,1) ;
var_sum=zeros(N-k,1);
% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y1(i)-y1bar)*(y2(i-k)-y2bar) * 1/(N-k);
   % cross_sum(i) = (y(i)-ybar)*(y(i-k)-ybar) * 1/(N); % original matlab
   %  implementation
end

% Denominator, unscaled variance
%yvarold = abs(y1-y1bar)'*abs(y2-y2bar)* 1/N ;

yvar1= (y1-y1bar)'*(y1-y1bar)* 1/N ; % variance y1
yvar2= (y2-y2bar)'*(y2-y2bar)* 1/N ; % variance y2

yvar=sqrt(yvar1*yvar2); % product of standard deviations

ta2 = sum(cross_sum) / yvar ;
% Numerator, unscaled covariance
% for i = (k+1):N
%     cross_sum(i) = (y1(i)-y1bar)*(y2(i-k)-y2bar) ;
%     var_sum(i)=abs(y1(i)-y1bar)*abs(y2(i-k)-y2bar);
%     
% end
% 
% 
% % divide 
% ta2 = sum(cross_sum) / sum(var_sum) ;


