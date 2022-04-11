function prs = setupDecode(gg,PriorHyp,win)
T = diff(win)+1;
ncells = length(gg);
dt = gg{1}.dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct inverse prior Hessian matrix and prior mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prs=mkbighess(PriorHyp,T);
prs.F = size(gg{1}.k,2);


% make time-interpolation matrix for downsamling spike-history current
MMntrp = makeInterpMatrix2(T,dt);


for i=1:ncells
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save history current and offset
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---  Compute spike times as integers in fine time sampling lattice -----
    r = ceil((gg{i}.tsp-dt*.001)/dt); % own spike train
    
    H = gg{i}.ihbas*gg{i}.ih;  % expand spike history filter
    finewin = ceil([win(1) win(2)+1-dt]./dt );
    yh = spikeconv_mex(r,H,finewin); % compute spike-history current
    yhc = MMntrp'*yh;%downsample spike-history current
    
      prs.bh{i} = yhc*dt + gg{i}.dc; % save history current and offset
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct STRF convolution matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    strf = flipud(gg{i}.k).';
    prs.K{i} = genkernhist2(T,T,[strf ]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save spike times as integers in bins equal to spectrogram bin size--
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rcourse = ceil(gg{i}.tsp).'; % own spike train
    rcourse=rcourse(win(1)<=rcourse & rcourse<=win(2)); % indices in correct window
    
    prs.spikes{i} = hist(rcourse,[win(1):win(2)]-0.5);
end

  