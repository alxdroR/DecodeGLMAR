function hhist = genkernhist2(T,N,h)
% hhist = genkernhist2(T,N,h)
% Given the  F X M spectro-temporal filter h, genkernhist2
% generates a T x N*F sparse stimulus history matrix.
% One can perform the following causal convolution with hhist
% 
% y_t = [hhist*x(:)]_t = sum_f sum_tau h(f,t-tau)*x(f,tau) 
%
% ADR

[F M] = size(h);
%hhist = zeros(T,N*F);
hhist = spalloc(T,N*F,1);

if T>=M
    % The first row of hhist will contain all temporal components of h
    t1star = T-M+1;
    if N >= t1star
        % loop through all blocks of h^T
        for tau=1:t1star
            hhist(tau:M+(tau-1),F*(tau-1)+1:F*tau) = h.';
        end
        if N >= M+t1star-1
            % loop through M
            for tau=1:M-1
                hhist(t1star+tau:T,F*(tau+t1star-1)+1:F*(tau+t1star)) = (h(:,1:M-tau)).';
            end
        else
            % loop through N*F
            for tau=1:N-t1star
                hhist(t1star+tau:T,F*(tau+t1star-1)+1:F*(tau+t1star)) = (h(:,1:M-tau)).';
            end
        end
    else
        % loop until N
        for tau=1:N
            hhist(tau:M+(tau-1),F*(tau-1)+1:F*tau) = h.';
        end
    end
else
    % The first row of hhist will NOT contain all temporal components of h
    if N >= T
        % loop through T
        for tau=1:T
            hhist(tau:T,F*(tau-1)+1:F*tau) = (h(:,1:T-tau+1)).';
        end
    else
        % loop through N*F
        for tau=1:N
            hhist(tau:T,F*(tau-1)+1:F*tau) = (h(:,1:T-tau+1)).';
        end
    end
end

 
 
