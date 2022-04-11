function [H , A] = spHessAR(ARcoeffs,T)
% [H , A] = spHessAR(ARcoeffs,T)
% spits out the sparse Hessian (Inverse covariance) of an AR process. 
% T is the duration of the time interval of interest  (size(H) = [T,T])
% ARcoeffs are the AR coeffs. 
% the Hessian is a toeplitz with nonzero (-m:m) diagonals (so it's (2m+1)-diagonal)
% where m is the order of the AR (m = length(ARcoeffs)). 
% A is the matrix such that A*x = ee; where x is the time series and 
% ee is a vector of white noises. 
% Yashar Ahmadian-2009

% X = toeplitz(zeros(ARorder,1),x);
% X = X';
% 
% k = (X'*X)\(X'x);
% X = toeplitz(zeros(ARorder,1),x);


ARcoeffs = reshape(ARcoeffs,[],1);
m = length(ARcoeffs);

if m<T
    A = sptoeplitz([1;-ARcoeffs;zeros(T-m-1,1)],[1,zeros(1,T-1)]);
else
    % no point here in making this sparse, really.
    A = sptoeplitz([1;-ARcoeffs(1:T-1)],[1,zeros(1,T-1)]);
end

H = A'*A;
