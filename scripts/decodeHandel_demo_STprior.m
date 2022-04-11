% Decode Gaussian spectrogram with spectrotemporal correlations 
% assuming a GLM likelihood and Gaussian prior.
%
% Code Blocks:  
%   1. Set up model params
%   2. Show simulated model responses to stimuli
%   3. Make training data (for fitting params to data) 
%   4. Fit GLM params via maximum-likelihood (requires optimization toolbox)
%      and fit prior hyperparameters (correlations in stimulus)
%      (demo default is to skip the costly ML step and just use true
%      parameters)
%   5. decode spectrogram
%   6. visualize results
 
clear all
global RefreshRate;  % Stimulus refresh rate (Stim frames per second)

%% 1.  Set parameters and display for GLM ============= % 

DTsim = .01; % Bin size for simulating model & computing likelihood.
nkt = 15;  % Number of time bins in temporal component of STRF;
ggsim = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
kt = ggsim.k*0;
kt(end-3) = 1;

% Make a spatial filter;
ncells = 150;  % Number of (non-coupled) cells to use in simulation.
nkx = 35; % number of frequency bins
xxk = [1:1:nkx]'; % index over frequency bins

kx = zeros(nkx,min(nkx,ncells));
for i=1:min(nkx,ncells)
    kx(:,i) = 1./sqrt(2*pi*4).*exp(-(xxk-i).^2/5); % construct spectral filter centered on bin i
end

[iht,ihbasOrthog,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,DTsim);
%ggsim.ih = ihbasis*[-10;-0.3;0;0;0]; % post-spike refractory period


figure(1);  % === Make Fig: model params =======================
Filt = kt*kx(:,1)'; % Make example space-time separable filter
ggsim.k = Filt./norm(Filt(:))*3; % Insert into simulation struct
ttk = [-nkt+1:0]';

subplot(3,5,[1,6]); % ------------------------------------------
plot(xxk,kx(:,1)); axis tight;
set(gca, 'xlim', [.5 nkx+.5]);
xlabel('space (pixels)');

subplot(3,5,[2,3,7,8]); % --------------------------------------
imagesc(ttk,xxk,ggsim.k');  axis xy
colormap gray;
axis off; 
title('example stimulus kernel k');

subplot(3,5,[12,13]); % ----------------------------------------
plot(ttk,kt);  axis tight;
%set(gca, 'xdir', 'reverse');
xlabel('time (frames)');



subplot(3,5,4:5); % --------------------------------------------
plot(ggsim.iht, ggsim.iht*0, 'k--', ggsim.iht, ggsim.ih); 
title('post-spike kernel h');
set(gca, 'xlim', ggsim.iht([1 end]),...
    'ylim',[min(ggsim.ih)*1.1 max(ggsim.ih)*1.5]);
subplot(3,5,9:10); % -------------------------------------------
plot(ggsim.iht, ihbasis);
axis tight;
xlabel('time after spike (frames)');
title('basis for h');

%% 2. Construct spectrogram & simulate the glm model response. ========= %

slen = 100; % Stimulus length (temporal bins)
swid = size(ggsim.k,2); % number of frequency bins

load('handel')
x = y(1e3:2e4);
window = 256;
nolap = round(max(0,(slen*window - length(x))/(slen-1)));
f = linspace(0,4e3,nkx); % frequency bins
[S,F,T,P] = spectrogram(x,window,nolap,f,Fs);


RefreshRate = 1/(T(2)-T(1));
Stim = 10*log10(P)';
Stim = (Stim-mean(Stim(:)))/(std(Stim(:)));
slen=length(T);

% w = hanning( (length(F)-1)*2,'periodic');
% overlap=4;
% [rec,D]=LSEE(S,w,overlap);

[tsp, vmem,Ispk] = simGLM(ggsim, Stim);

% ==== Make Figure ========
figure(2); 
tt = [DTsim:DTsim:slen]'/RefreshRate;
subplot(221); %------------------------
imagesc(T,F,Stim'); axis xy
title('Handel Normalized PSD');
xlabel('time (s)');
ylabel('frequency(Hz)');

subplot(223); %------------------------
plot(tt, vmem-Ispk, tt, Ispk, 'r', T(ceil(tsp)), .1*ones(size(tsp)), 'r.');
title('stim- and spike-induced currents'); axis tight;
xlabel('time (s)');
subplot(222); %------------------------
plot(tt, vmem, T(ceil(tsp)), max(vmem)*ones(size(tsp)), 'ro');
title('net voltage (dots = spike times)'); 
axis tight;
% -----Run repeat simulations ------
nrpts = 3;        % number of repeats to draw
subplot(224);
for j = 1:nrpts;
    [tsp1,vmem1] = simGLM(ggsim,Stim);
    plot(tt,vmem1, 'color', rand(3,1)); hold on;
end
axis tight; hold off;
title('repeat responses to same stim');
xlabel('time (s)');
set(gcf,'Name','Example Cell Simulation Results')

sound(x,Fs)
%% 3. Generate some training data ========================================
slen = 300;
window = 256;
nolap = round(max(0,(slen*window - length(x))/(slen-1)));
f = linspace(0,4e3,nkx); % frequency bins
[S,F,T,P] = spectrogram(x,window,nolap,f,Fs);


RefreshRate = 1/(T(2)-T(1));
Stim = 10*log10(P)';
Stim = (Stim-mean(Stim(:)))/(std(Stim(:)));
slen=length(T);


% sample ncells  
tsp = cell(ncells,1);
nsp = zeros(ncells,1);
kt = zeros(nkt,ncells);
for i=1:ncells
    
    % index for assigning best frequency to cell i
    fidx = mod(i,nkx);
    if fidx == 0
        fidx = nkx;
    end
   % tk = [0:nkt-1]';
   % b1 = rand*0.3+0.3; b2 = 2*b1;
   % k1 = 1/(gamma(6)*b1)*(tk/b1).^5 .* exp(-tk/b1);  % Gamma pdfn
   % k2 = 1/(gamma(6)*b2)*(tk/b2).^5 .* exp(-tk/b2);  % Gamma pdf
  %  kt(:,i) = flipud(k1-k2./1.5);
    kt(end-5,i)=1;
    Filt = kt(:,i)*kx(:,fidx)'; % Make space-time separable filter
    ggsim.k = Filt./norm(Filt(:))*5; % Insert into simulation struct
    [tsp{i} vk{i}] = simGLM(ggsim,Stim);  % run model
    nsp(i) = length(tsp{i});
end


%% estimate spectrogram posterior parameters: GLM coefficients, spectrotemporal correlations of song spectrogram
 
%  Fit GLM (traditional version) via max likelihood

%  FOR DEMONSTRATION PURPOSES FITTING IS TURNED OFF 
%  AND WE WILL DECODE ASSUMING OUR FITS PRODUCE TRUE 
%  PARAMETERS

dofitting = 0; % If this boolean variable is false -- assumes GLM parameters are known when and set to true values

gg1 = cell(ncells,1);
negloglival = zeros(ncells,1);
if ~dofitting
    dc= ggsim.dc;
    ihc = ihbasOrthog\ggsim.ih;
end
    

for i=1:ncells
    if dofitting
        % Compute STA and use as initial guess for k
        sta0 = simpleSTC(Stim,tsp{1},nkt);
        sta = reshape(sta0,nkt,[]);
    else
        sta = zeros(nkt,swid);
    end
    %  Initialize params for fitting --------------
    gg0 = makeFittingStruct_GLM(sta,DTsim);
    gg0.tsp = tsp{i};
    gg0.tspi = 1;
    
    if dofitting
        opts = {'display', 'iter', 'maxiter', 100};
        % Do ML estimation of model params
        [gg1{i}, negloglival(i)] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)
    else
        gg1{i}=gg0;
        gg1{i}.dc = dc;
        gg1{i}.ih = ihc;
        fidx = mod(i,nkx);
        if fidx == 0
            fidx = nkx;
        end
        Filt = kt(:,i)*kx(:,fidx)'; % Make space-time separable filter
        gg1{i}.k = Filt./norm(Filt(:))*3; % Insert into simulation struct
        
    end
end

%  Determine prior parameters -- temporal correlations parameterized by 
%  autoregressive process of order m and spectral correlations empircally fit

m = 26;
PriorHyp=HP_SC(Stim',0,m);
%% decode spectrogram

% decode spectrogram from bin --  win(1) to win(3)
T = 233;
win = [4 4-1+T];

% OPTIMIZATION PARAMETERS
options=optimset('GradObj','on', 'Hessian','on','Display','on','MaxIter',400,'TolX',1e-9);


% set up parameters to evaluate stimulus posterior
prs=setupDecode(gg1,PriorHyp,win);



% Handle for posterior function with current input parameters
pos = @(x) logpos_stim_matrix(prs,x);


x0 = prs.priormean;
% call Newton-Raphson optimization routine
[x fn grad Hdmb ef] = NewtRaph(pos,x0,options);



%% visualize decoded results 
figure(3);  % === Training data =======================
SP=zeros(max(nsp),ncells);

for i=1:ncells
    SP(1:length(tsp{i}),i) = tsp{i};
end
Tmax=200; nmax = 100;
subplot(311); % ------------------------------------------
Stimseg = Stim(win(1):win(2),:)';
mx = max(abs(Stimseg(:)));
imagesc(Stimseg,[-mx mx]); axis xy
axis off; 
title('stimulus');

subplot(312); % ----------------------------------------
plotraster(SP(:,1:min(nmax,ncells)));
xlabel('Temporal bin number')
ylabel('Cell number')
title('Raster Plot')

subplot(313)
xresh = reshape(x,swid,T);
imagesc(xresh,[-mx mx]) ; axis xy
axis off; 
title('decoded stimulus');


