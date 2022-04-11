function [lpneg g H pg] = logpos_stim_matrix(prs,x)
% [lpneg g H pg] = logpos_stim_matrix(prs,x)
% Evaluates stimulus dependent part of GLM log posterior with a Gaussian prior at input.
% Also returns gradient and Hessian at input point.
% 
% INPUT:
% x  - is a D X 1 stimulus vector.  To construct this vector from the 
% spectrogram that gives the power at time t and frequency f, x(f,t),
% use x(:)
%
% prs - structure containing GLM and prior parameters
% 
%       prs.spikes{k}(i,j) stores a 1(0) when cell k fires/(doesn't fire) a spike at bin i on
%       trial j.
% 
%       prs.dt - temporal size of a spike bin (in seconds).
% 
%       prs.K{j} - STRF convolution matrix for cell j. Use genkernhist2(T,T,strf)
%       [prs.K{j}*x(:)]_t = sum_f sum_tau prs.K{j}(f,t-tau)*x(f,tau) 
%
%       prs.bh{k}(i,t) is a matrix that gives cell k's bias plus
%       cell k's spike-history filter convolved with the spike-train from trial i
%       at time t.
%
%       prs.priormean - D x 1 vector storing  prior mean
% 
%       prs.hessx - D x D prior covariance matrix.  (usually sparse)
%
% OUTPUT:
% lpneg - THe negative of the stimulus dependent part of the log posterior 
% evaluated at input value x.
%
% g - The negative gradient with respect to (w.r.t) the stimulus of the log
% posterior evaluated at the input x
%
% H - The negative Hessian matrix w.r.t the stimulus of the log posterior
% evaluated at the input
% 
% pg - prior gradient evaluated at the input.
%


% adr 
% 11-30-09

global RefreshRate

dt = 1/RefreshRate;
% Total dimension of quantity being decoded
D = length(x);
% find the number of frequencies the STRF filters
nf = prs.F;
% number of time-bins
T = D/nf;

% total number of different cell types
ncells = length(prs.bh);


lpneg = 0;
if nargout > 1
    g = zeros(D,1);
end
if nargout > 2
    H = spalloc(D,D,1);
end
% loop over cells
for j=1:ncells
    
    K = prs.K{j};
    ystim = K*x;
    
    
    niid = size(prs.spikes{j},1);
    % loop over i.i.d cells
    for i =1:niid
        % bias and history terms
        Sbh = prs.bh{j};
        % intensity fnc
        s = Sbh + ystim;
        f = exp(s)*dt;
        
        %calculate the negative of the log posterior
        lpneg = lpneg - (-sum(f) + prs.spikes{j}(i,:)*s);
        
        %if number of outputs is greater than 1, calculate gradient
        if nargout>1
            g = g + (-1)*K.'*(-f + prs.spikes{j}(i,:).' );
        end
        
        %also, if number of outputs is greater than 2, calculate hessian
 
        if nargout>2
            HL =spdiags(-f,0,T,T); 
            Hj = -1*K.'*HL*K;
             H = H + Hj;
        end
    end

end

%% add prior contribution
lprior = -0.5*(x-prs.priormean).'*prs.hessx*(x-prs.priormean);
lpneg = lpneg -lprior;

if nargout>1
    pg =prs.hessx*(x-prs.priormean); 
   % g = g + prs.hessx*(x-prs.priormean);
     g = g + pg;
end

if nargout>2
    H = H +  prs.hessx;
end

end




