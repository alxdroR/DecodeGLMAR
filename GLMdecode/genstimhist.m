function X = genstimhist(M,x)
% X = genstimhist(M,x)
% INPUT: 
% x - double - F x N stimulus matrix 
% % M is the maximum lagtime of the convolution
% 
% OUTPUT:
% X - double - N X M*F stimulus history matrix
% 
% Given the F x N stimulus matrix x, genstimhist
% generates an N x M*F stimulus history matrix, X.
% X is constructed so that 
% y = Xh
% is the convolution of h with x.

[F N] = size(x);
X = zeros(N,M*F);
if N<M
    for tau =1:N
        % legal values in convolution sum j
        j0 = max(1,tau+1-M);
        jmax = min(tau,N);
        j = jmax:-1:j0;
        X(tau,1:length(j)) = x(j);
    end
else
    for p=0:M-1
        X(p+1:end,F*p+1:F*(p+1)) = x(:,1:end-p).';
    end
end
