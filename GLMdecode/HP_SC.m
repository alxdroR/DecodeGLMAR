function songHyp = HP_SC(X,whiten,M)
% songHyp = HP_SC(X,whiten,M)
%
%  Fit  spectro-temporal correlations to spectrogram in X assuming 
%  seperability.  Fit an AR process of order M to temporal components 
%  and empirically determine frequency matrix.
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Samples we will use for parameter estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = size(X,2);
mus = mean(X,2);
Y = X - mus*ones(1,T);
% Fit frequency statistics
Cfwhiten = cov(Y.');
CfInv = Cfwhiten\eye(size(Cfwhiten));



if whiten
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Whiten stimulus: The following whitenting procedure allows to treat
    % the temporal variations at each frequency as independent realizations
    % from a process that only depends on time, not frequency.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform a whitening transformation on the spectrogram
    
    Y = sqrtm(CfInv)*Y;
end

% Fit temporal statistics
AST = Y';
[am err] = arburg(AST(:),M);
b = sqrt(err);

if ~whiten
    % Since our model assumes a separable covariance matrix
    % we need to scale the matched variance.  Currently err
    % is such that TEMPORAL variance proportional to variance in song,
    %
    % err = factor*var(AST(:));
    %
    % but we want it to match variance of
    %AS2 = Y'*sqrtm(CfInv)';
    %tempvar = var(AS2(:));
    tempvar = trace(Cfwhiten)/length(Cfwhiten);
    
    % so b will be set to
    b = sqrt(tempvar*err/(var(AST(:))));
end

% prior Hyperparameters
songHyp.a = am;
songHyp.b = b;
songHyp.ssvar = b;
songHyp.freqinv = CfInv;
songHyp.mean =  mus;


