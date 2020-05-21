function [L,P,f,Lc] = spod(X,varargin)
%SPOD Spectral proper orthogonal decomposition
%   [L,P,F] = SPOD(X) returns the spectral proper orthogonal decomposition
%   of the data matrix X whose first dimension is time. X can have any
%   number of additional spatial dimensions or variable indices.
%   The columns of L contain the modal energy spectra. P contains the SPOD
%   modes whose spatial dimensions are identical to those  
%   of X. The first index of P is the frequency and
%   the last one the mode number ranked in descending order 
%   by modal energy. F is the frequency vector. If DT is not specified, a
%   unit frequency sampling is assumed. For real-valued data, one-sided spectra
%   are returned. Although SPOD(X) automatically chooses default spectral 
%   estimation parameters, it is strongly encouraged to manually specify
%   problem-dependent parameters on a case-to-case basis.
%
%   [L,P,F] = SPOD(X,WINDOW) uses a temporal window. If WINDOW is a vector, X
%   is divided into segments of the same length as WINDOW. Each segment is
%   then weighted (pointwise multiplied) by WINDOW. If WINDOW is a scalar,
%   a Hamming window of length WINDOW is used. If WINDOW is omitted or set
%   as empty, a Hamming window is used.
%
%   [L,P,F] = SPOD(X,WINDOW,WEIGHT) uses a spatial inner product weight in
%   which the SPOD modes are optimally ranked and orthogonal at each
%   frequency. WEIGHT must have the same spatial dimensions as X. 
%
%   [L,P,F] = SPOD(X,WINDOW,WEIGHT,NOVERLAP) increases the number of
%   segments by overlapping consecutive blocks by NOVERLAP snapshots.
%   NOVERLAP defaults to 50% of the length of WINDOW if not specified. 
%
%   [L,P,F] = SPOD(X,WINDOW,WEIGHT,NOVERLAP,DT) uses the time step DT
%   between consecutive snapshots to determine a physical frequency F. 
%
%   [L,P,F] = SPOD(XFUN,...,OPTS) accepts a function handle XFUN that 
%   provides the i-th snapshot as x(i) = XFUN(i). Like the data matrix X, x(i) can 
%   have any dimension. It is recommended to specify the total number of 
%   snaphots in OPTS.nt (see below). If not specified, OPTS.nt defaults to 10000.
%   OPTS.isreal should be specified if a two-sided spectrum is desired even
%   though the data is real-valued, or if the data is initially real-valued, 
%   but complex-valued-valued for later snaphots.
%
%   [L,P,F] = SPOD(X,WINDOW,WEIGHT,NOVERLAP,DT,OPTS) specifies options:
%   OPTS.savefft: save FFT blocks to avoid storing all data in memory [{false} | true]
%   OPTS.deletefft: delete FFT blocks after calculation is completed [{true} | false]
%   OPTS.savedir: directory where FFT blocks and results are saved [ string | {'results'}]
%   OPTS.savefreqs: save results for specified frequencies only [ vector | {all} ]
%   OPTS.loadfft: load previously saved FFT blocks instead of recalculating [{false} | true]
%   OPTS.mean: provide a mean that is subtracted from each snapshot [ array of size X | 'blockwise' | {temporal mean of X; 0 if XFUN} ]
%   OPTS.nsave: number of most energtic modes to be saved [ integer | {all} ]
%   OPTS.isreal: complex-valuedity of X or represented by XFUN [{determined from X or first snapshot if XFUN is used} | logical ]
%   OPTS.nt: number of snapshots [ integer | {determined from X; defaults to 10000 if XFUN is used}]
%   OPTS.conflvl: confidence interval level [ scalar between 0 and 1 | {0.95} ]
%   OPTS.normvar: normalize each block by pointwise variance [{false} | true]
%
%   [L,PFUN,F] = SPOD(...,OPTS) returns a function PFUN instead of the SPOD
%   data matrix P if OPTS.savefft is true. The function returns the j-th
%   most energetic SPOD mode at the i-th frequency as p = PFUN(i,j) by
%   reading the modes from the saved files. Saving the data on the hard
%   drive avoids memory problems when P is large. The modal energy spectra
%   are returned in L, and saved in a separate file 'spod_energy.mat'.
%
%   [L,P,F,Lc] = SPOD(...) returns the confidence interval Lc of L. By
%   default, the lower and upper 95% confidence levels of the j-th
%   most energetic SPOD mode at the i-th frequency are returned in
%   Lc(i,j,1) and Lc(i,j,2), respectively. The OPTS.conflvl*100% confidence
%   interval is returned if OPTS.conflvl is set. For example, by setting 
%   OPTS.conflvl = 0.99 we obtain the 99% confidence interval. A 
%   chi-squared distribution is used, i.e. we assume a standard normal 
%   distribution of the SPOD eigenvalues.
%
%   References:
%     [1] Towne, A. and Schmidt, O. T. and Colonius, T., Spectral proper 
%         orthogonal decomposition and its relationship to dynamic mode
%         decomposition and resolvent analysis, J. of Fluid Mech. 847, 821–867, 2018
%     [2] Schmidt, O. T. and Towne, A. and Rigas, G. and Colonius, T. and 
%         Bres, G. A., Spectral analysis of jet turbulence, J. of Fluid Mech. 855, 953–982, 2018
%     [3] Lumley, J. L., Stochastic tools in turbulence, Academic Press, 
%         1970
%
% O. T. Schmidt (oschmidt@ucsd.edu), A. Towne, T. Colonius
% Last revision: 20-May-2020

if nargin==6
    opts = varargin{5};
    if ~isfield(opts,'savefft')
        opts.savefft = false;
    end
else
    opts.savefft = false;
end

% get problem dimensions
if isa(X, 'function_handle')
    xfun    = true;
    if isfield(opts,'nt')
        nt  = opts.nt;
    else
        warning(['Please specify number of snapshots in "opts.nt".' ...
        ' Trying to use default value of 10000 snapshots.'])
        nt  = 10000;
    end
    sizex   = size(X(1));
    nx      = prod(sizex);
    dim     = [nt sizex];
else 
    xfun    = false;
    dim     = size(X);
    nt      = dim(1);
    nx      = prod(dim(2:end));
end

% Determine whether data is real-valued or complex-valued-valued to decide on one- or two-sided
% spectrum. If "opts.isreal" is not set, determine from data. If data is
% provided through a function handle XFUN and opts.isreal is not specified,
% deternime complex-valuedity from first snapshot.
if isfield(opts,'isreal')
    isrealx = opts.isreal;
elseif ~xfun
    isrealx = isreal(X);
else
    x1      = X(1);
    isrealx = isreal(x1);
    clear x1;
end

% get default spectral estimation parameters and options
[window,weight,nOvlp,dt,nDFT,nBlks] = spod_parser(nt,nx,isrealx,varargin{:});

% determine correction for FFT window gain
winWeight   = 1/mean(window);

% Use data mean if not provided through "opts.mean".
blk_mean    = false;
if isfield(opts,'mean')
    if strcmp('blockwise',opts.mean)
        blk_mean    = true;       
    end
end
switch xfun
    case true
        if isfield(opts,'mean')&&(~blk_mean)
            x_mean      = opts.mean(:);
            mean_name   = 'user specified';
        elseif blk_mean
            mean_name   = 'blockwise mean';
        else
            x_mean      = 0;
            warning('No mean subtracted. Consider providing long-time mean through "opts.mean" for better accuracy at low frequencies.');            
            mean_name   = '0';
        end
    case false
        if blk_mean
            mean_name   = 'blockwise mean';
        else
            x_mean      = mean(X,1);
            x_mean      = x_mean(:);
            mean_name   = 'long-time (true) mean';
        end
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

% set default for confidence interval 
if nargout==4
    confint = true;
    if ~isfield(opts,'conflvl')
        opts.conflvl = 0.95;
    end
    xi2_upper   = 2*gammaincinv(1-opts.conflvl,nBlks);
    xi2_lower   = 2*gammaincinv(  opts.conflvl,nBlks);
    Lc          = zeros(nFreq,nBlks,2);
else
    confint = false;
    Lc      = [];
end

% set defaults for options for saving FFT data blocks
if ~isfield(opts,'savefreqs'),  opts.savefreqs  = 1:nFreq;   end
if opts.savefft
    if ~isfield(opts,'savedir'),    opts.savedir   = pwd;        end
    saveDir = fullfile(opts.savedir,['nfft' num2str(nDFT) '_novlp' num2str(nOvlp) '_nblks' num2str(nBlks)]);
    if ~exist(saveDir,'dir'),       mkdir(saveDir);              end
    if ~isfield(opts,'nsave'),      opts.nsave      = nBlks;     end   
    if ~isfield(opts,'deletefft'),  opts.deletefft  = true;      end    
    omitFreqIdx = 1:nFreq; omitFreqIdx(opts.savefreqs) = []; 
end 
if ~isfield(opts,'loadfft'),    opts.loadfft    = false;     end
if ~isfield(opts,'normvar'),    opts.normvar    = false;     end

% loop over number of blocks and generate Fourier realizations
disp(' ')
disp('Calculating temporal DFT')
disp('------------------------------------')
if ~opts.savefft
    Q_hat = zeros(nFreq,nx,nBlks);
end
Q_blk       = zeros(nDFT,nx);
for iBlk    = 1:nBlks
    
    % check if all FFT files are pre-saved
    all_exist   = 0;
    if opts.loadfft        
        for iFreq = opts.savefreqs
            if ~isempty(dir(fullfile(saveDir,['fft_block' num2str([iBlk iFreq],'%.4i_freq%.4i.mat')])))
                all_exist   = all_exist + 1;
            end
        end
    end
    if opts.loadfft && all_exist==length(opts.savefreqs)
        disp(['loading FFT of block ' num2str(iBlk) '/' num2str(nBlks) ' from file'])
    else
        
        % get time index for present block
        offset                  = min((iBlk-1)*(nDFT-nOvlp)+nDFT,nt)-nDFT;
        timeIdx                 = (1:nDFT) + offset;
        disp(['block ' num2str(iBlk) '/' num2str(nBlks) ' (' ...
            num2str(timeIdx(1)) ':' num2str(timeIdx(end)) ')'])
        % build present block
        if blk_mean, x_mean = 0; end
        if xfun
            for ti = timeIdx
                x                  = X(ti);
                Q_blk(ti-offset,:) = x(:) - x_mean;
            end
        else
            Q_blk   = bsxfun(@minus,X(timeIdx,:),x_mean.');
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
        
        % window and Fourier transform block
        Q_blk                   = bsxfun(@times,Q_blk,window);
        Q_blk_hat               = winWeight/nDFT*fft(Q_blk);
        Q_blk_hat               = Q_blk_hat(1:nFreq,:);
        % correct Fourier coefficients for one-sided spectrum
        if isrealx
            Q_blk_hat(2:end-1,:)    = 2*Q_blk_hat(2:end-1,:);
        end
        if ~opts.savefft
            % keep FFT blocks in memory
            Q_hat(:,:,iBlk)         = Q_blk_hat;
        else
            for iFreq = opts.savefreqs
                file = fullfile(saveDir,['fft_block' num2str([iBlk iFreq],'%.4i_freq%.4i')]);
                Q_blk_hat_fi        = single(Q_blk_hat(iFreq,:));
                save(file,'Q_blk_hat_fi','-v7.3');
            end
        end
    end
end

% loop over all frequencies and calculate SPOD
L   = zeros(nFreq,nBlks);
disp(' ')
disp('Calculating SPOD')
disp('------------------------------------')
if ~opts.savefft
    % keep everything in memory (default)
    P   = zeros(nFreq,nx,nBlks);
    for iFreq = 1:nFreq
        disp(['frequency ' num2str(iFreq) '/' num2str(nFreq) ' (f=' num2str(f(iFreq),'%.3g') ')'])
        Q_hat_f             = squeeze(Q_hat(iFreq,:,:));
        M                   = Q_hat_f'*bsxfun(@times,Q_hat_f,weight)/nBlks;
        [Theta,Lambda]      = eig(M);
        Lambda              = diag(Lambda);
        [Lambda,idx]        = sort(Lambda,'descend');
        Theta               = Theta(:,idx);
        Psi                 = Q_hat_f*Theta*diag(1./sqrt(Lambda)/sqrt(nBlks));
        % Theta               = Theta*diag(sqrt(Lambda)*sqrt(nBlks));
        P(iFreq,:,:)        = Psi;
        L(iFreq,:)          = abs(Lambda);
        if confint
            Lc(iFreq,:,1)       = L(iFreq,:)*2*nBlks/xi2_lower;
            Lc(iFreq,:,2)       = L(iFreq,:)*2*nBlks/xi2_upper;
        end
    end
    P   = reshape(P,[nFreq dim(2:end) nBlks]);
else
    % save FFT blocks on hard drive (for large data)
    for iFreq = opts.savefreqs
        disp(['frequency ' num2str(iFreq) '/' num2str(nFreq) ' (f=' num2str(f(iFreq),'%.3g') ')'])
        % load FFT data from previously saved file
        Q_hat_f             = zeros(nx,nBlks);
        for iBlk    = 1:nBlks
            file    = fullfile(saveDir,['fft_block' num2str([iBlk iFreq],'%.4i_freq%.4i')]);
            load(file,'Q_blk_hat_fi');
            Q_hat_f(:,iBlk) = Q_blk_hat_fi;
        end
        M                   = Q_hat_f'*bsxfun(@times,Q_hat_f,weight)/nBlks;
        [Theta,Lambda]      = eig(M);
        Lambda              = diag(Lambda);
        [Lambda,idx]        = sort(Lambda,'descend');
        Theta               = Theta(:,idx);
        Psi                 = Q_hat_f*Theta*diag(1./sqrt(Lambda)/sqrt(nBlks));
        % Theta               = Theta*diag(sqrt(Lambda)*sqrt(nBlks));
        if opts.nsave>0
            Psi             = single(reshape(Psi(:,1:opts.nsave),[dim(2:end) opts.nsave]));
            file            = fullfile(saveDir,['spod_f' num2str(iFreq,'%.4i')]);
            save(file,'Psi','-v7.3');
        else
            Psi             = [];
        end
        L(iFreq,:)          = abs(Lambda);
        if confint
            Lc(iFreq,:,1)       = L(iFreq,:)*2*nBlks/xi2_lower;
            Lc(iFreq,:,2)       = L(iFreq,:)*2*nBlks/xi2_upper;          
        end
    end
    % return anonymous function that loads SPOD modes from files instead of
    % actual solution
    P   = @(iFreq,iMode) getmode(iFreq,iMode,saveDir);
    
    % return sparse matrices
    L   = sparse(L);
    Lc  = sparse(Lc);
    
    file            = fullfile(saveDir,'spod_energy');
    save(file,'L','Lc','f','-v7.3');
    
    % clean up
    if opts.deletefft
        for iFreq = opts.savefreqs
            for iBlk    = 1:nBlks
                file    = fullfile(saveDir,['fft_block' num2str([iBlk iFreq],'%.4i_freq%.4i.mat')]);
                delete(file);
            end
        end
    end
    
    disp(' ')
    disp(['Results saved in folder ' saveDir])
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = getmode(iFreq,iMode,saveDir)
%GETMODE Anonymous function that loads SPOD modes from files
file        = matfile(fullfile(saveDir,['spod_f' num2str(iFreq,'%.4i')]));
dim         = size(file.Psi);
for di = 1:length(dim)-1
    idx{di} = 1:dim(di);
end
idx{di+1}   = iMode;
P           = file.Psi(idx{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window,weight,nOvlp,dt,nDFT,nBlks] = spod_parser(nt,nx,isrealx,varargin)
%SPOD_PARSER Parser for SPOD parameters

% read input arguments from cell array
window = []; weight = []; nOvlp = []; dt = [];
nvarargin = length(varargin);

if nvarargin >= 1
    window = varargin{1};
    if nvarargin >= 2
        weight   = varargin{2};
        if nvarargin >= 3
            nOvlp   = varargin{3};
            if nvarargin >= 4
                dt      = varargin{4};
            end
        end
    end
end

window = window(:); weight = weight(:);

% check arguments and determine default spectral estimation parameters
% window size and type
if isempty(window)
    nDFT        = 2^floor(log2(nt/10));
    window      = hammwin(nDFT);
    window_name = 'Hamming';
elseif length(window)==1
    nDFT        = window;
    window      = hammwin(window);
    window_name = 'Hamming';
elseif length(window) == 2^nextpow2(length(window))
    nDFT        = length(window);
    window_name = 'user specified';
else
    nDFT        = length(window);
    window_name = 'user specified';
end

% block overlap
if isempty(nOvlp)
    nOvlp = floor(nDFT/2);
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
if nDFT < 4 || nBlks < 2
    error('Spectral estimation parameters not meaningful.');
end
    
% display parameter summary
disp(' ')
disp('SPOD parameters')
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window] = hammwin(N)
%HAMMWIN Standard Hamming window of lenght N
    window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
end
