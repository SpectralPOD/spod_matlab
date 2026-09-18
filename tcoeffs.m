function [a] = tcoeffs(X,P,window,weight,nModes,St_ind,useParallel)
%TCOEFFS Continuously-discrete temporal expansion coefficients of SPOD modes
%   A = TCOEFFS(X,P,WINDOW,WEIGHT,NMODES) returns the
%   continuously-discrete temporal SPOD mode expansion coefficients of the
%   leading NMODES modes. P is the data matrix of SPOD modes returned by
%   SPOD. X, WINDOW and WEIGHT are the same variables as for SPOD. If
%   WINDOW is a scalar, a Hamming window of length WINDOW is used. If
%   WEIGHT is empty, a uniform weighting of 1 is used.
%
%   A = TCOEFFS(X,P,WINDOW,WEIGHT,NMODES,ST_IND) restricts the computation
%   to the inclusive frequency-index range ST_IND = [FIRST LAST]. An empty
%   ST_IND selects all frequencies. The first dimension of A follows the
%   selected frequency range.
%
%   A = TCOEFFS(X,P,WINDOW,WEIGHT,NMODES,ST_IND,USEPARALLEL) uses a
%   PARFOR loop when USEPARALLEL is true. Parallel processing is disabled
%   by default so that the standard call has no toolbox dependencies.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency-time analysis, low-rank
%         reconstruction and denoising of turbulent flows using SPOD,
%         Journal of Fluid Mechanics 926, A26, 2021
%
% A. Nekkanti (aknekkan@eng.ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
% Revision history:
%   7-Oct-2022: Brandon Yeung <byeung@ucsd.edu>
%   26-May-2023: Tianyi Chu <tic173@eng.ucsd.edu> reordered the
%                computation to project spatially before the temporal FFT.
%   13-Aug-2026: Joel Weightman <joel.weightman@monash.edu> added block-wise
%                parallel processing, vectorized frequency processing, and
%                support for selecting a frequency-index range.
%   20-Aug-2026: Integrated the direct-projection update while preserving
%                subsequent real/complex mode-array handling.
%   3-Sep-2026:  Made parallel execution explicitly opt-in.

dims        = size(X);
nt          = dims(1);
nGrid       = prod(dims(2:end));
isRealX     = isreal(X);

window      = window(:);
weight      = weight(:);

% default window size and type
if isscalar(window)
    window  = hammwin(window);
end

nDFT        = length(window);
nHalf       = ceil(nDFT/2);
winCorr_fac = 1/(mean(window)*nDFT);

% inner product weight
if isempty(weight)
    weight  = ones(nGrid,1);
elseif ~isscalar(weight) && numel(weight)~=nGrid
    error('tcoeffs:InvalidWeight', ...
        'WEIGHT must be empty, scalar, or contain one value per grid point.')
end

X           = reshape(X,nt,nGrid);
meanX       = mean(X,1);

if isRealX
    nFreq   = ceil(nDFT/2)+1;
else
    nFreq   = nDFT;
end

if size(P,1)~=nFreq
    error('tcoeffs:InvalidModes', ...
        'The first dimension of P must contain %d frequencies.',nFreq)
end
if nModes<1 || nModes~=fix(nModes) || nModes>size(P,ndims(P))
    error('tcoeffs:InvalidModeCount', ...
        'NMODES must be an integer between 1 and the number of modes in P.')
end

if nargin<6 || isempty(St_ind)
    St_ind  = [1 nFreq];
elseif numel(St_ind)~=2 || any(St_ind~=fix(St_ind)) || ...
        St_ind(1)<1 || St_ind(2)>nFreq || St_ind(1)>St_ind(2)
    error('tcoeffs:InvalidFrequencyRange', ...
        'ST_IND must be an inclusive range [FIRST LAST] within 1:%d.',nFreq)
end
freqInd     = St_ind(1):St_ind(2);
nFreqOut    = numel(freqInd);

if nargin<7 || isempty(useParallel)
    useParallel = false;
elseif ~isscalar(useParallel) || ...
        (~islogical(useParallel) && ~ismember(useParallel,[0 1]))
    error('tcoeffs:InvalidParallelOption', ...
        'USEPARALLEL must be a scalar logical value.')
end

P_dims      = size(P);
P           = permute(P,[1 length(P_dims) 2:length(P_dims)-1]);
P           = reshape(P,nFreq,P_dims(end),nGrid);
P           = P(freqInd,1:nModes,:);
P_proj      = reshape(permute(P,[3 2 1]),nGrid,nModes*nFreqOut);

% Linear indices select FFT bin FREQIND(k) from the projection onto the
% modes at frequency FREQIND(k), for every selected frequency and mode.
fftInd      = (1:nModes).' + (0:nFreqOut-1)*nModes;
fftInd      = fftInd + (freqInd-1)*nModes*nFreqOut;
window      = reshape(window,1,1,nDFT);

% Preserve the original double-precision output contract. Input data are
% projected in blocks so that the low-dimensional intermediate remains
% bounded for long time series.
a           = complex(zeros(nFreqOut,nModes,nt));
maxBlockSize= 5e4;
nBlocks     = ceil(nt/maxBlockSize);

% Correction for windowing and zero-padding. These are based on global
% time indices, so block boundaries do not affect the result.
corr        = ones(nt,1);
windowVec   = window(:);
for i=1:nt
    if i<nHalf+1
        corr(i) = 1/(winCorr_fac*sum(windowVec(nHalf-i+1:nDFT)));
    elseif i>nt-nHalf+1
        corr(i) = 1/(winCorr_fac*sum(windowVec(1:nt+nHalf-i)));
    end
end

disp(' ')
disp('Calculating expansion coefficients')
disp('------------------------------------')

for iBlock=1:nBlocks
    iFirst      = (iBlock-1)*maxBlockSize+1;
    iLast       = min(iBlock*maxBlockSize,nt);
    ntBlock     = iLast-iFirst+1;

    % Each output uses NDFT consecutive samples from the globally padded
    % signal. Only the real samples needed by this block are projected.
    dataFirst   = max(1,iFirst-nHalf);
    dataLast    = min(nt,iLast-nHalf+nDFT-1);
    insertFirst = dataFirst-(iFirst-nHalf)+1;
    insertLast  = insertFirst+dataLast-dataFirst;

    fprintf('block %d/%d (snapshots %d:%d)\n', ...
        iBlock,nBlocks,iFirst,iLast)

    X_block     = X(dataFirst:dataLast,:)-meanX;
    X_proj_data = P_proj'*(weight.*transpose(X_block));
    X_proj_data = reshape(X_proj_data,nModes,nFreqOut,[]);
    X_proj      = zeros(nModes,nFreqOut,ntBlock+nDFT-1, ...
                        'like',X_proj_data);
    X_proj(:,:,insertFirst:insertLast) = X_proj_data;

    aBlock      = complex(zeros(nFreqOut,nModes,ntBlock));
    corrBlock   = corr(iFirst:iLast);
    if useParallel
        parfor i=1:ntBlock
            a_fft         = fft(X_proj(:,:,i:i+nDFT-1).*window,nDFT,3);
            aBlock(:,:,i) = corrBlock(i)*winCorr_fac* ...
                reshape(a_fft(fftInd),nModes,nFreqOut).';
        end
    else
        for i=1:ntBlock
            a_fft         = fft(X_proj(:,:,i:i+nDFT-1).*window,nDFT,3);
            aBlock(:,:,i) = corrBlock(i)*winCorr_fac* ...
                reshape(a_fft(fftInd),nModes,nFreqOut).';
        end
    end

    a(:,:,iFirst:iLast) = aBlock;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window] = hammwin(N)
%HAMMWIN Standard Hamming window of length N
    window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
end
