function [prs Ctinv Cfinv Ct] = mkbighess(PriorHyp,T)
% [prs Ctinv Cfinv Ct] = mkbighess(PriorHyp,m,T)
% Construct input parameters for a Gaussian distribution with a covariance
% matrix composed of the tensor product of a user input matrix and the covariance 
% matrix of an autoregressive process with parameters and duration specified by user.
% Gaussian covariance is normalized so that variance is equal to the
% variance of the user input covariance matrix.
% 
% Output is stored in a format compatible with scripts for computing 
% GLM + Gaussian log-posterior
% 
% INPUT:
% PriorHyp - structure containing input matrix and autoregressive parameters
%         The process, z(n), obeys the equation
%         PriorHyp.b.order{m}z(n) =  PriorHyp.a.order{m}(1)z(n-1) +
%         PriorHyp.a.order{m}(2)z(n-2) + ...  PriorHyp.a.order{m}(m)z(n-m)
%
%         PriorHyp.freqinv - inverse of FxF user input matrix.
% 
%         PriorHyp.mean - F dimensional user input vector specifying
%         time-independent Gaussian mean.
%
% m - order of Autoregressive process
% 
% T - total number of time bins for AR process.
%
% OUTPUT:
% prs - output structure containing Gaussian parameters
%       prs.hessx - Inverse of (TF x TF) Gaussian covariance matrix
%       prs. priormean - TF dimensional vector specifying time-independent 
%       Gaussian mean
%       prs.F - number of dimensions of user specified input matrix.
% 
% Ctinv - Inverse of AR process covariance matrix
% 
% Cfinv  - Inverse of user specified input matrix
% 
% Ct -  AR covariance matrix
%
% Calls the following functions
% See also: spHessAR, kron

% ADR
% updated 03/11/2012 

burn = 200;
Ctinv = spHessAR(-1*PriorHyp.a(2:end),T+burn);
% ensure that steady-state variance of process matches song-variance
Ctinv2 = Ctinv/(PriorHyp.b^2);

% Remove burn-time and start at steady-state
Ct = Ctinv2\eye(T+burn);
Ct = Ct(burn+1:end,burn+1:end);
Ctinv3 = Ct\eye(T);
Ctinv3(abs(Ctinv3) < 1e-9) = 0;
Ctinv = sparse(Ctinv3);


Cfinv = PriorHyp.freqinv*(trace(Ct)/T);
prs.hessx = kron(Ctinv,Cfinv); % Kron will maintain the sparsity inherent in Ctinv


mu = PriorHyp.mean*ones(1,T);
prs.priormean = mu(:);

end
