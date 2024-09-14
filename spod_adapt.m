function [L,P,f,A,nWin] = spod_adapt(X,varargin)
%SPOD_ADAPT Adaptive sine-taper SPOD by Yeung & Schmidt[1]
%   [L,P,F] = SPOD_ADAPT(X) returns the spectral proper orthogonal decomposition
%   of the data matrix X whose first dimension is time. X can have any
%   number of additional spatial dimensions or variable indices. The
%   columns of L contain the modal energy spectra. P contains the SPOD
%   modes whose spatial dimensions are identical to those of X. The first
%   index of P is the frequency and the last one the mode number ranked in
%   descending order by modal energy. F is the frequency vector. If DT is
%   not specified, a unit frequency sampling is assumed. For real-valued
%   data, adjusted one-sided eigenvalue spectra are returned. Although
%   SPOD(X) automatically chooses default spectral estimation parameters,
%   the user is encouraged to manually specify problem-dependent parameters
%   on a case-to-case basis.
%
%   [L,P,F] = SPOD_ADAPT(X,WINDOW) uses a temporal window. If WINDOW is a
%   vector, X is divided into segments of the same length as WINDOW. Each
%   element in WINDOW is the number of sine windows to be used at the
%   corresponding frequency. If WINDOW is a scalar, one sine window of
%   length WINDOW is used. If WINDOW is omitted or set as empty, three sine
%   windows of length equal to the time dimension of X is used.
%   Multitaper-Welch estimates are computed when the number of windows is
%   greater than one; see [4,5] for details. Standard SPOD is computed otherwise.
%
%   [L,P,F] = SPOD_ADAPT(X,WINDOW,WEIGHT) uses a spatial inner product weight in
%   which the SPOD modes are optimally ranked and orthogonal at each
%   frequency. WEIGHT must have the same spatial dimensions as X. 
%
%   [L,P,F] = SPOD_ADAPT(X,WINDOW,WEIGHT,NOVERLAP) increases the number of
%   segments by overlapping consecutive blocks by NOVERLAP snapshots.
%   NOVERLAP defaults to 0 if not specified. Multitaper SPOD neither
%   requires nor benefits from overlap.
%
%   [L,P,F] = SPOD_ADAPT(X,WINDOW,WEIGHT,NOVERLAP,DT) uses the time step DT
%   between consecutive snapshots to determine a physical frequency F. 
%
%   [L,P,F] = SPOD_ADAPT(X,WINDOW,WEIGHT,NOVERLAP,DT,OPTS) specifies options:
%   OPTS.savefreqs: store results for specified frequencies only [ vector | {all} ]
%   OPTS.mean: provide a mean that is subtracted from each snapshot [ array of size X | 'blockwise' | {temporal mean of X} ]
%   OPTS.nsave: number of most energtic modes to be stored [ integer | {all} ]
%   OPTS.isreal: complex-valuedity of X [{determined from X} | logical ]
%   OPTS.normvar: normalize each block by pointwise variance [{false} | true]
%   OPTS.compress: lossless time-domain compression for compute and memory efficiency [{false} | true]
%   OPTS.truncate: number of ranks in compressed data to keep [integer | {all}]
%   OPTS.adaptive: determine number of windows using adaptive algorithm [{false} | true]
%   OPTS.tol: convergence tolerance for adaptive algorithm [integer | {1e-6}]
%
%   [L,P,F,A] = SPOD_ADAPT(...) returns the block-wise expansion coefficients in A.
%
%   [L,P,F,A,nWin] = SPOD_ADAPT(...) returns the window numbers in nWin.
%
%   Reference:
%     [1] Yeung, B. C. Y., Schmidt, O. T. Adaptive spectral proper orthogonal 
%         decomposition of broadband-tonal flows, Theor. Comput. Fluid Dyn. 38, 355–374, 2024,
%         DOI 10.1007/s00162-024-00695-0
%
% B. Yeung (byeung@ucsd.edu) and O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 10-Sep-2024

if nargin==6
    opts = varargin{5};
else
    opts = [];
end

% get problem dimensions
dim     = size(X);
nt      = dim(1);
nx      = prod(dim(2:end));

% check whether data is single-precision
if isa(X,'single'), precision = 'single'; else, precision = 'double'; end

% Determine whether data is real-valued or complex-valued to decide on one- or two-sided
% spectrum. If "opts.isreal" is not set, determine from data. If data is
% provided through a function handle XFUN and opts.isreal is not specified,
% deternime complex-valuedity from first snapshot.
if isfield(opts,'isreal')
    isrealx = opts.isreal;
else
    isrealx = isreal(X);
end

% get default spectral estimation parameters and options
[weight,nOvlp,dt,nDFT,nBlks,nWin] = spod_parser(nt,nx,isrealx,varargin{:});

% Use data mean if not provided through "opts.mean".
blk_mean    = false;
if isfield(opts,'mean')
    if strcmp('blockwise',opts.mean)
        blk_mean    = true;       
    end
end

if blk_mean
    mean_name   = 'blockwise mean';
elseif isfield(opts,'mean')
    x_mean      = opts.mean(:);
    mean_name   = 'user specified';
else    
    x_mean      = mean(X,1);
    x_mean      = x_mean(:);
    mean_name   = 'data mean';
end
disp(['Mean                      : ' mean_name]);

% obtain frequency axis
f = (0:nDFT-1)/dt/nDFT;
if isrealx
    f = (0:ceil(nDFT/2))/nDFT/dt;
else
    if mod(nDFT,2)==0
        f(nDFT/2+1:end) = f(nDFT/2+1:end)-1/dt;
    else
        f((nDFT+1)/2+1:end) = f((nDFT+1)/2+1:end)-1/dt;
    end
end
nFreq = length(f);

% sine window central difference parameters
wrapIdx = @(idx,N) 1+mod(idx-1,N);
sinetapNorm = sqrt(2/(nDFT+1));

% set defaults for options
if ~isfield(opts,'savefreqs'),              opts.savefreqs  = 1:nFreq;   end
if ~isfield(opts,'normvar'),                opts.normvar    = false;     end
if ~isfield(opts,'compress'),               opts.compress   = false;     end
if ~isfield(opts,'truncate'),               opts.truncate   = 0;         end
if isempty(nWin) && ~isfield(opts,'tol'),   opts.tol        = 1e-6;      end

% switch off compression for small data
if opts.compress && nx<=nt
    opts.compress = false;
    warning('Spatial dimensions smaller than time dimension. Ignoring compression.')
end

% time-domain compression
if opts.compress
    X = X(:,:).' - x_mean;
    sqrtW = sqrt(weight);
    X = sqrtW.*X;   
    if (nx < 8*nt) && ~opts.truncate
        % if spatial dimension is small and truncation is unneeded, compress using QR
        [U,X] = qr(X,0); nx = size(U,2);
    else
        % if spatial dimension is large or truncation is needed, compress using EVD
        M = X'*X;
        nx = size(M,1);
        [Theta,Lambda]      = eig(M,'vector');
        [Lambda,idx] = sort(Lambda,'descend');
        Theta = Theta(:,idx);
        Sigma = sqrt(Lambda);
        U = X; clear X;
        U = U*(Theta*diag(1./Sigma));
        X = diag(Sigma)*Theta';
    end
    X = X.';
end

% loop over number of blocks and generate Fourier realizations
disp(' ')
disp('Calculating temporal DFT')
disp('------------------------------------')
Q_hat = zeros(nDFT*2,nx,nBlks,precision);
for iBlk    = 1:nBlks
        
        % get time index for present block
        offset                  = min((iBlk-1)*(nDFT-nOvlp)+nDFT,nt)-nDFT;
        timeIdx                 = (1:nDFT) + offset;
        % build present block
        if blk_mean, x_mean = 0; end
        if ~opts.compress
            Q_blk   = bsxfun(@minus,X(timeIdx,:),x_mean.');
        else
            Q_blk   = X(timeIdx,:);
        end
        % if block mean is to be subtracted, do it now that all data is
        % collected
        if blk_mean
            Q_blk   = bsxfun(@minus,Q_blk,mean(Q_blk,1));           
        end        
        
        % normalize by pointwise variance
        if opts.normvar
            Q_var   = sum(bsxfun(@minus,Q_blk,mean(Q_blk,1)).^2,1)/(nDFT-1);
            % address division-by-0 problem with NaNs             
            Q_var(Q_var<4*eps)  = 1; 
            Q_blk   = bsxfun(@rdivide,Q_blk,Q_var);
        end        
        
        disp(['block ' num2str(iBlk) '/' num2str(nBlks) ' (snapshots ' ...
            num2str(timeIdx(1)) ':' num2str(timeIdx(end)) ')'])
        % Fourier transform block using zero-padded FFT
        Q_blk = cat(1,Q_blk,zeros(nDFT,nx,precision));
        Q_blk = sqrt(dt)*fft(Q_blk);
        Q_hat(:,:,iBlk)         = Q_blk;
end
clear Q_blk Q_blk_hat

% losslessly compressed data can optionally be truncated, resulting in lossy compression
Q_hat = permute(Q_hat,[2 3 1]);
if opts.compress && opts.truncate
    Q_hat_full = Q_hat;
    nx = opts.truncate;
    Q_hat = Q_hat(1:nx,:,:);
end

% get adaptive windows
if isempty(nWin)
    nWin = adaptiveWindows(Q_hat,wrapIdx,sinetapNorm,opts);
end

% always compute final results on untruncated data
if opts.compress && opts.truncate
    Q_hat = Q_hat_full;
    nx    = size(U,2);
end

% loop over all frequencies and calculate SPOD
disp(' ')
if max(nWin)>1
    disp('Computing multitaper SPOD')
else
    disp('Computing SPOD')
end
disp('------------------------------------')
nSamples    = nBlks*max(nWin);
if ~isfield(opts,'nsave'), opts.nsave = nBlks*min(nWin(1:opts.savefreqs(end))); end
L   = zeros(nFreq,nSamples,precision);
A   = zeros(nFreq,nSamples,nSamples,precision);
P   = zeros(opts.savefreqs(end),nx,opts.nsave,precision);
for iFreq = opts.savefreqs
    nTapers = nWin(iFreq);
    nSamples = nBlks*nTapers;
    
    % parabolic weights penalize high-order sine windows
    taperWeights = parabWeights(nTapers);
    taperWelchWeights = repmat(taperWeights.',nBlks,1)/nBlks;
    taperWelchWeights = taperWelchWeights(:).';
    
    disp(['frequency ' num2str(iFreq) '/' num2str(nFreq) ' (f=' num2str(f(iFreq),'%.3g') ') with ' num2str(nTapers) ' windows'])
    
    % recover windowed DFT from unwindowed DFT using sine window central difference
    Q_blk_hat_fi_l  = Q_hat(:,:,wrapIdx((iFreq-1)*2+1-(1:nTapers),nDFT*2));
    Q_blk_hat_fi_u  = Q_hat(:,:,wrapIdx((iFreq-1)*2+1+(1:nTapers),nDFT*2));
    Q_hat_f         = sinetapNorm*(Q_blk_hat_fi_l-Q_blk_hat_fi_u)/(2i);
    
    Q_hat_f = reshape(Q_hat_f,[nx nSamples]);    
    Q_hat_f = Q_hat_f.*sqrt(taperWelchWeights);
    
    if ~opts.compress
        M                   = Q_hat_f'*bsxfun(@times,Q_hat_f,weight);
    else
        M                   = Q_hat_f'*Q_hat_f;
    end
    [Theta,Lambda]      = eig(M);
    Lambda              = diag(Lambda);
    [Lambda,idx]        = sort(Lambda,'descend');
    Theta               = Theta(:,idx);
    Psi                 = Q_hat_f*Theta(:,1:opts.nsave)*diag(1./sqrt(Lambda(1:opts.nsave))/sqrt(nSamples));
    P(iFreq,:,:)        = Psi(:,1:opts.nsave);
    A(iFreq,1:nSamples,1:nSamples)        = diag(sqrt(Lambda))*Theta'.*sqrt(1./taperWelchWeights);
    L(iFreq,1:nSamples)          = abs(Lambda);
    % correct for one-sided spectrum
    if isrealx && iFreq~=1 && iFreq~=nFreq
        L(iFreq,:) = 2*L(iFreq,:);
    end
end
clear Q_hat

% recover SPOD modes from compression basis
if opts.compress
    nx = prod(dim(2:end));
    P = permute(P,[2 3 1]);
    P = U*P(:,:)./sqrtW;
    P = reshape(P,[nx opts.nsave opts.savefreqs(end)]);
    P = permute(P,[3 1 2]);
end
P   = reshape(P,[opts.savefreqs(end) dim(2:end) opts.nsave]);

% return the same number of eigenvalues and coefficients for all frequencies
A   = A(:,1:nBlks*min(nWin(1:opts.savefreqs(end))),:);
L   = L(:,1:nBlks*min(nWin(1:opts.savefreqs(end))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [weight,nOvlp,dt,nDFT,nBlks,nWin] = spod_parser(nt,nx,isrealx,varargin)
%SPOD_PARSER Parser for SPOD parameters

% read input arguments from cell array
window = []; weight = []; nOvlp = []; dt = []; adaptive = false;
nvarargin = length(varargin);

if nvarargin >= 1
    window = varargin{1};
    if nvarargin >= 2
        weight   = varargin{2};
        if nvarargin >= 3
            nOvlp   = varargin{3};
            if nvarargin >= 4
                dt      = varargin{4};
                if nvarargin >= 5
                    opts    = varargin{5};
                    if isfield(opts,'adaptive')
                        adaptive = opts.adaptive;
                    end
                end
            end
        end
    end
end

% check window contains only integers
if ~isempty(window)
    assert(sum(window~=round(window))==0,'Window must be one of: integer scalar | integer vector | [].')
end

% check arguments and determine default spectral estimation parameters
% window size and type
window_name     = 'sine';
if isempty(window)
    nDFT        = nt;
    nWin        = 3*ones(nDFT,1);
elseif length(window)==1
    nDFT        = window;
    nWin        = ones(nDFT,1);
elseif length(window)>=2
    nDFT        = length(window);
    nWin        = window;
end

if adaptive
    nWin        = [];
end

weight      = weight(:);

% block overlap
if isempty(nOvlp)
    nOvlp = 0;
elseif nOvlp > nDFT-1
    error('Overlap too large.')
end

% time step between consecutive snapshots
if isempty(dt)
    dt = 1;
end

% inner product weight
if isempty(weight)
    weight      = ones(nx,1);
    weight_name = 'uniform';
elseif numel(weight) ~= nx
    error('Weights must have the same spatial dimensions as data.');
else
    weight_name = 'user specified';
end

% number of blocks
nBlks = floor((nt-nOvlp)/(nDFT-nOvlp));

% test feasibility
if ~isempty(nWin)
    if nDFT < 4 || nBlks*min(nWin(nWin>0)) < 3
        error('Spectral estimation parameters not meaningful.');
    end
end
    
% display parameter summary
disp(' ')
if ~isempty(nWin)
    if max(nWin)==1
        disp('SPOD parameters')
    else
        disp('Multitaper SPOD parameters')
    end
else
    disp('Adaptive SPOD parameters')
end
disp('------------------------------------')
if isrealx
    disp('Spectrum type             : one-sided (real-valued signal)')
else
    disp('Spectrum type             : two-sided (complex-valued signal)')
end
disp(['No. of snaphots per block : ' num2str(nDFT)])
disp(['Block overlap             : ' num2str(nOvlp)])
disp(['No. of blocks             : ' num2str(nBlks)])
disp(['Windowing fct. (time)     : ' window_name])
disp(['Weighting fct. (space)    : ' weight_name])
if adaptive
    disp(num2str(opts.tol,'No. of windows            : adaptive tol = %g'))
elseif max(nWin)==min(nWin)
    disp(['No. of windows            : ' num2str(nWin(1))])
else
    disp(num2str([min(nWin) max(nWin)],'No. of windows            : varies between %i and %i'))
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu_k] = parabWeights(K)
%PARABWEIGHTS Parabolic weights corresponding to K sine tapers
    Nk = K*(4*K-1)*(K+1)/6;
    k = (1:K).';
    mu_k = 1/Nk*(K^2-(k-1).^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TF = anynan(A)
%ANYNAN checks if A contains NaNs
TF = sum(isnan(A))>0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function similarity = calcModeSimilar(q1,q2)
%MODESIMILAR calculates the similarity between complex SPOD modes q1 and q2

q1 = double(q1(:));
q2 = double(q2(:));

similarity = dot(q1,q2)/norm(q1)/norm(q2);
similarity = abs(similarity);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nWin = adaptiveWindows(Q_hat,wrapIdx,sinetapNorm,opts)
%ADAPTIVEWINDOWS implements the adaptive algorithm and returns a frequency-dependent nWin

% get data parameters
nx      = size(Q_hat,1);
nBlks   = size(Q_hat,2);
nDFT    = size(Q_hat,3)/2;

disp(' ')
disp('Adaptively converging nWin')
disp('------------------------------------')

% store only the leading mode to track mode convergence
P = zeros(opts.savefreqs(end),nx,'like',Q_hat);

% initialize window numbers
nWin            = nan(size(opts.savefreqs));
allConverged    = false;
if nBlks>=3, nTapers_init = 1; else, nTapers_init = 3; end
nTapers         = nTapers_init - 1;

% increment window numbers until all frequencies are converged to within tolerance
while ~allConverged
    nTapers = nTapers+1;
    if nTapers>nTapers_init
        P_prev = P;
    end
    unconvergedFreqs = opts.savefreqs(isnan(nWin));
    disp(['Converging ' num2str(length(unconvergedFreqs)) '/' num2str(length(opts.savefreqs)) ' frequencies using ' num2str(nTapers) ' windows'])
    
    % SPOD for unconverged frequencies
    for iFreq = unconvergedFreqs
        nSamples = nBlks*nTapers;
        taperWeights = parabWeights(nTapers);
        taperWelchWeights = repmat(taperWeights.',nBlks,1)/nBlks;
        taperWelchWeights = taperWelchWeights(:).';
        
        Q_blk_hat_fi_l = Q_hat(:,:,wrapIdx((iFreq-1)*2+1-(1:nTapers),nDFT*2));
        Q_blk_hat_fi_u = Q_hat(:,:,wrapIdx((iFreq-1)*2+1+(1:nTapers),nDFT*2));
        Q_hat_f   = sinetapNorm*(Q_blk_hat_fi_l-Q_blk_hat_fi_u)/(2i);
        Q_hat_f = reshape(Q_hat_f,[nx nSamples]);        
        Q_hat_f             = Q_hat_f.*sqrt(taperWelchWeights);
        
        if nSamples<=nx
            M                   = Q_hat_f'*Q_hat_f;
        else
            M                   = Q_hat_f*Q_hat_f';
        end
        
        if size(M,1)<60
            [Theta,Lambda]      = eig(M);
            Lambda              = diag(Lambda);
            [Lambda,idx]        = max(Lambda);
            Theta               = Theta(:,idx);
        else
            [Theta,Lambda]      = eigs(double(M),1);
        end
        if nSamples<=nx
            Psi                 = Q_hat_f*Theta*diag(1./sqrt(Lambda)/sqrt(nSamples));
        else
            Psi                 = Theta/sqrt(nSamples);
        end
        
        P(iFreq,:)        = Psi;
        
        % check if mode is converged
        if nTapers>nTapers_init
            similarity = calcModeSimilar(P_prev(iFreq,:),P(iFreq,:));
            if (1-similarity)<=opts.tol
                nWin(iFreq) = nTapers;
            else
                if iFreq<opts.savefreqs(end)
                    if nTapers-nWin(iFreq+1)==1
                        nWin(iFreq) = nTapers;
                    end
                end
                if iFreq>1
                    if nTapers-nWin(iFreq-1)==1
                        nWin(iFreq) = nTapers;
                    end
                end
            end
        end
    end
    if ~anynan(nWin)
        allConverged = true;
    end
end
end