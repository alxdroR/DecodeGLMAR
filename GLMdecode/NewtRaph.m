function [x fn grad THTA A eflag] = NewtRaph(f,x0,varargin)
% function [x fn grad A] = NewtRaph(f,x0,options)
% implements the conjugate-gradient algorithm for minimizing a quadratic
% function, f, in N variables given the Hessian function, H.  
%
% INPUT: 
%   f - Handle to function being evaluated. IMPORTANT: F MUST 
%   RETURN FUNCTION VALUE, GRADIENT, AND HESSIAN IN THAT ORDER
%   x0 - starting point (column vector)
%   options - structure of options
% OUTPUT:
%   x - minimum point
%   fn - value of function at minimum
%   grad - gradient of function at minimum
%   A - Hessian of function at minimum
%
% Using the notation from Numerical Recipies, Second Edition(pg 516)

% ADR
% ?

g = x0;
h = x0;
 THTA(:,1) = g;
% check maximum number of iterations
if isempty(varargin)
    options.MaxIter = 800;
    options.TolFun = 1e-8;
    options.Display = 'on';
else
    options = varargin{1};
end

if ~isempty(options.MaxIter)
    maxit = options.MaxIter;
else
    maxit = 100;
end

options.TolFun = [];
options.Display = 'on';


% Function tolerance
if ~isempty(options.TolFun)
    tol = options.TolFun;
else
    tol = 1e-8;
end
gradtol = 1e-6;
maxnsteps = 70;
trouble = 0;
% evaluate function at starting point
[f0 grad A] = f(g);
%f0 = -f0;
%grad = -grad;
%A = -A;
fnstuck = 0;
fn(1) = f0;
if ~isempty(options.Display) &  strcmp(options.Display,'on')
    verbose = 1;
else
    verbose = 0;
end
%    fprintf('fval %d  \n',fn(1));
if verbose
fprintf('%10s %10s %15s  %15s %15s\n','Iteration','FunEvals','Function Val','OptCond','CPU time');
end

use1Dsolver =0 ;
if use1Dsolver
    opts1D = optimset('Gradobj','on','Hessian','on','display','final','DerivativeCheck','off','maxiter',1,'Largescale','on');
end

for i=1:maxit
    tic;
    % initial point
    g0 = g;
    
    %c = cond(A);
    %if c > 1e4
     %   'Large Condition number'
    %    [U D] = eigs(A);
        %gdir=pinv(A,mean(real(diag(-D))))*grad;
   % else
        
    % direction we move
    
    gdir = A\grad;
    %end
    alpha = 1;
    % should be a while loop?
    
    % current value of function
    ftest0 = fn(i);
    
    if ~use1Dsolver
        for stepsize =1:maxnsteps
            % step
            g = g0 - alpha*gdir;
            ftest = f(g);
            % break when we move to a spot where fnc is lower
            if ftest < ftest0
                eflag = 1;
                break
            else
                alpha = alpha*0.5;
            end
        end
        
        
%         if stepsize == maxnsteps
%        %     fprintf('1D solver didnt budge in %d steps: \n ftest0 = %d, \n ftest(alpha=%d) = %d',maxnsteps,ftest0,alpha,ftest)
%             trouble = 1;
%         end
    else
        %astar = fminunc(@(alpha) loss1D(alpha,g0,gdir),1,opts1D);
        [f2 g2 H2] = loss1D(alpha,g0,gdir);
        astar = alpha - g2/H2;
        g = g0-astar*gdir;
    end
    % update A
    [fn(i+1) grad A] = f(g);
    %fn(i+1) = -fn(i+1);
    fn(i+1) = fn(i+1);
    %   grad = -grad;
    %  A = -A;
    %if ( abs(fn(i+1) - fn(i)) < 1e-5 )
    %    g(end)
   % end
    % check if change in function is below tolerance
    opcon = sum(abs(grad));
    if ( abs(fn(i+1) - fn(i)) < tol || opcon < eps )
        if norm(grad) > 9e-5
        %    fprintf('\n fval difference = %d, but ||g|| = %d\n....this could take awhile \n',abs(fn(i+1) - fn(i)),norm(grad));
            
            % Record the stuck value of 
            fnstuck0 = fn(i+1);
         %   if(trouble || abs(fnstuck0 - fnstuck) < tol)
                fprintf('\nProgram has failed to satisfy all criteria:\n Function change = %d but \n ||g|| = %d \n',abs(fn(i+1) - fn(i)),norm(grad));
                x = g;
                 THTA(:,i+1) = g;
                eflag = 3;
                return;
          %  end
          
        else
             cput=toc;
           fprintf('%10d %10d %15.5e %15.5e %15.5e\n',i,stepsize,fn(i+1),opcon,cput);
            THTA(:,i+1) = g;
            x=g;
            eflag = 3;
            return;
        end
    elseif ~isempty(options.Display) &  strcmp(options.Display,'on')
        cput=toc;
        fprintf('%10d %10d %15.5e %15.5e %15.5e\n',i,stepsize,fn(i+1),opcon,cput);
        
        THTA(:,i+1) = g;
    end
end
x = g;
eflag =0;
warning('Maximum number of function evaluations reached');

end

function [neglogli dL H] = loss1D(x,t0,tdir)
% [neglogli, dL, H] = loss1D(x,t0,tdir)
%
% Compute negative log-likelihood of data under the GLM model
% for the special case where the input parameter vector, theta,
% is of the form
% 
% theta = t0 + x*tdir
% 
% and t0,tdir are known
% Uses arbitrary nonlinearity 'nlfun' instead of exponential
%
% Inputs:
%   t0 = varies depending on which components are held fixed
%         Without any fixed components 
%        = [kprs - weights for stimulus kernel
%          dc   - dc current injection
%          ihprs - weights on post-spike current];
%
%          To fix a component, just leave the appropriate 
%          t0 elements blank.
% Outputs:
%      neglogli = negative log likelihood of spike train
%      dL = gradient with respect to x
%      H = hessian

etol = 1e-100;

% global vars for optimization
global MSTM MMntrp SPNDS SPNDS2 OPRS RefreshRate indEXT

% Extract some vals from OPRS (Opt Prs);
nprs = 1; % 
nt = OPRS.nt;
nx = OPRS.nx; 
nktot = nx*nt;
nh = OPRS.nh; 
nh2 = OPRS.nh2;
nsp = length(SPNDS);
nCoupled = length(SPNDS2);
nhtot = nh+nh2*nCoupled;
dt = OPRS.dt; ndt = round(1/dt);
slen = size(MSTM,1); rlen = slen*ndt;

prs{1} = t0;prs{2}=tdir;
% Unpack GLM prs;
for j=1:2
kprs{j} = prs{j}(indEXT.kt);
dc{j} = prs{j}(indEXT.dc);
ihprs{j} = prs{j}([indEXT.ih indEXT.ih2]);
end
% determine who is held fixed
%fix = find([isempty(indEXT.kt(:)),isempty(indEXT.dc),isempty(indEXT.ih(:)),isempty(indEXT.ih2(:)) ]);

% Initialize likelihood, gradient and Hessian -----------
neglogli = 0;
dL = zeros(nprs,1);
H = zeros(nprs,nprs);

for jch = 1:OPRS.nchunks
    %------- Compute indices relevant to this chunk ---------------
    iwin = OPRS.ichunk([jch, jch+1]); % abs fine bin win
    iwin1 = [iwin(1)+1,iwin(2)];
    jwin = OPRS.jchunk(jch,:);   % abs window of stimulus, coarse
    jj = [1,diff(jwin)];
    ilen = diff(iwin);  
    jlen = diff(jwin);
    ii = [(iwin(1)-jwin(1)*ndt+1),(iwin(1)-jwin(1)*ndt+ilen)];
    %   if (iwin(1)/ndt-1)  == jwin(1)  %--- evenly-dividing bins ----
    %   ii = [ndt+1, ndt+ilen];
    
    for j = 1:2
        % -------- Compute stim filter reponse -----------------------
        
        if ~isempty(kprs{j})
            SS = MSTM(jwin(1)+1:jwin(2),:);  % Relevant stim
            MM = MMntrp(ii(1):ii(2),jj(1):jj(2)); % Interpolation filter for this chunk
            ystm = SS*kprs{j};  % filtered stimulus
            ystmhi(:,j) = MM*ystm;
            if ~isempty(dc{j})
                ystmhi(:,j) = ystmhi(:,j) + dc{j};
            end
        else
            if ~isempty(dc{j})
                ystmhi(:,j) = dc{j}*ones(ilen,1);
            else
                ystmhi(:,j) = zeros(ilen,1);
            end
        end
    end
    
    % -------------- Compute net h current -----------------------
    
    if OPRS.ihflag
        Ih = zeros(ilen,nh+nh2*nCoupled);
        Ih(:,1:nh) = spikeconv_mex(SPNDS, OPRS.ihbas, iwin1);
        for jcpl = 1:nCoupled
            Ih(:,nh+nh2*(jcpl-1)+1:nh+nh2*jcpl) = ...
                spikeconv_mex(SPNDS2{jcpl}, OPRS.ihbas2, iwin1);
        end
        for j=1:2
            yh(:,j) = Ih*ihprs{j};
            Iinj(:,j) = ystmhi(:,j)+yh(:,j);
        end
    else
        for j=1:2
            Iinj(:,j) = ystmhi(:,j);
        end
    end
    
    % -------- Extract fixed linear component -----------------------
    Iinj(:,1) = Iinj(:,1) + OPRS.fixedyi{jch};
    
    
    % ---------  Compute likelihood itself  ------------------------
    i_sp = in(SPNDS, iwin+.5)-iwin(1);
    [rr,drr,ddrr] = OPRS.nlfun(Iinj(:,1) + x*Iinj(:,2));
    spflag = ~isempty(i_sp);

    % Check for zero-values of rr
    iiz = find(rr <= etol);
    rr(iiz) = etol; % Set value to small
    drr(iiz) = 0;  ddrr(iiz) = 0;  % Set derivs here to 0

    Trm1 = sum(rr)*dt/RefreshRate;  % non-spike term
    if spflag, 
        Trm2 = -sum(log(rr(i_sp))); % spike term
    else, Trm2 = 0;
    end    
    neglogli = neglogli + (Trm1 + Trm2);

    % ---------  Compute Gradient -----------------
    if (nargout > 1)
        
        % Non-spiking terms
        dL0 = drr'*Iinj(:,2);
        
         if spflag
            frac1 = drr(i_sp)./rr(i_sp);
            dL1 = frac1'*Iinj(i_sp,2);
        else
            dL1 = 0;
        end
        dLdx = dL0*dt/RefreshRate - dL1;
        dL = dL + dLdx;
        
    end

    % ---------  Compute Hessian -----------------
    if nargout > 2
        H0 = ddrr'*(Iinj(:,2).^2)*dt/RefreshRate;
       
        if spflag  % -------  If spikes in window ------
            frac2 = (rr(i_sp).*ddrr(i_sp) - drr(i_sp).^2)./rr(i_sp).^2;
            H1 = frac2'*(Iinj(i_sp,2).^2);
            Hx = H0 - H1;
        end
        H = H + Hx;
    end
    
end
end