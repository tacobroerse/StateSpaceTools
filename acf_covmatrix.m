function ta = acf_covmatrix(C,p)
% ACF - Compute Autocorrelations Through p Lags when covariance matrix is
% known
% >> myacf = acf_covmatrix(C,p) 
%
% Inputs:
% C - covariance matrix, nxn
% p - total number of lags, 1x1 integer
%
% Output:
% myacf - px1 vector containing autocorrelations
%        (First lag computed is lag 1. Lag 0 not computed)
%
%
% A bar graph of the autocorrelations is also produced, with
% rejection region bands for testing individual autocorrelations = 0.
%
% Note that lag 0 autocorelation is not computed, 
% and is not shown on this graph.
%
% Example:
% >> acf_covmatrix(randn(100,1), 10)
%


% --------------------------
% USER INPUT CHECKS
% --------------------------

[n1, n2] = size(C) ;
if n2 ~= n1
    error('Input matrix C must be an nxn column vector')
end

[a1, a2] = size(p) ;
if ~((a1==1 & a2==1) & (p<n1))
    error('Input number of lags p must be a 1x1 scalar, and must be less than size of matrix C')
end



% -------------
% BEGIN CODE
% -------------

ta = zeros(p,1) ;
global N 
N = max(size(C)) ;
%global ybar 
%ybar = mean(C); 

% Collect ACFs at each lag i
for i = 1:p
   ta(i) = acf_k(C,i) ; 
end

% % Plot ACF
% % Plot rejection region lines for test of individual autocorrelations
% % H_0: rho(tau) = 0 at alpha=.05
% bar(ta)
% line([0 p+.5], (2.32)*(1/sqrt(N))*ones(1,2))
% line([0 p+.5], (-2.32)*(1/sqrt(N))*ones(1,2))
% 
% % Some figure properties
% line_hi = (2.32)*(1/sqrt(N))+.05;
% line_lo = -(2.32)*(1/sqrt(N))-.05;
% bar_hi = max(ta)+.1 ;
% bar_lo = -max(ta)-.1 ;
% 
% if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
%     axis([0 p+.60 line_lo line_hi])
% else
%     axis([0 p+.60 bar_lo bar_hi])
% end
% title({' ','Sample Autocorrelations',' '})
% xlabel('Lag Length')
% set(gca,'YTick',[-1:.20:1])
% % set number of lag labels shown
% if (p<28 & p>4)
%     set(gca,'XTick',floor(linspace(1,p,4)))
% elseif (p>=28)
%     set(gca,'XTick',floor(linspace(1,p,8)))
% end
% set(gca,'TickLength',[0 0])




% ---------------
% SUB FUNCTION
% ---------------
function ta2 = acf_k(C,k)
% ACF_K - Autocorrelation at Lag k
% acf(C,k)
%
% Inputs:
% C - covariance matrix
% k - which lag to compute acf
% 

global N
cki = zeros(N-k,1) ;

% Numerator, covariance
for i = 1:N-k
    cki(i) = C(i,i+k);

end

% summation and normalisation
ck = 1/(N-k)*sum(cki);

% Denominator, variance
c0 = 1/N*sum(diag(C));

ta2 = ck/c0 ;

