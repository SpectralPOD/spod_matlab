function [a] = tcoeffs(X,P,window,weight,nModes)
%TCOEFFS Continuously-discrete temporal expansion coefficients of SPOD
%modes
%   [A] = TCOEFFS(X,P,WINDOW,WEIGHT,NMODES) returns the
%   continuously-discrete temporal SPOD mode expansion coefficients of the
%   leading NMODES modes. P is the data matrix of SPOD modes returned by
%   SPOD. X, WINDOW and WEIGHT are the same variables as for SPOD. If no
%   windowing function is specified, a Hamming window of length WINDOW will
%   be used. If WEIGHT is empty, a uniform weighting of 1 is used.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency–time analysis, low-rank 
%         reconstruction and denoising of turbulent flows using SPOD, 
%         Journal of Fluid Mechanics 926, A26, 2021
%
% A. Nekkanti (aknekkan@eng.ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 7-Oct-2022 Brandon Yeung <byeung@ucsd.edu>

dims        = size(X);
nt          = dims(1);
nGrid       = prod(dims(2:end));

window = window(:); weight = weight(:);

% default window size and type
if length(window)==1
    window  = hammwin(window);
end

nDFT        = length(window);
winWeight   = 1/mean(window);

% inner product weight
if isempty(weight)
    weight  = ones(nGrid,1);
end

X           = reshape(X,nt,nGrid);
X           = X-mean(X,1);              % subtract mean

ndims       = size(P);
P           = permute(P,[1 length(ndims) 2:length(ndims)-1]);
if isreal(X)
    nFreq   = ceil(nDFT/2)+1;
    P       = reshape(P,ceil(nDFT/2)+1,ndims(end),nGrid);
else
    nFreq   = nDFT;
    P       = reshape(P,nDFT,ndims(end),nGrid);
end

weight     = reshape(weight,1,nGrid);
% zero-padding
X          = [zeros(ceil(nDFT/2),nGrid); X; zeros(ceil(nDFT/2),nGrid);];
a          = zeros(nFreq,nModes,nt);
winCorr_fac= winWeight/nDFT;
disp(' ')
disp('Calculating expansion coefficients')
disp('------------------------------------')
for i=1:nt
    X_blk               = fft(X(i:i+nDFT-1,:).*window);
    X_blk               = X_blk(1:nFreq,:);
    % correction for windowing and zero-padding
    if (i<ceil(nDFT/2)+1)
        corr 	= 1/(winCorr_fac*sum(window(ceil(nDFT/2)-i+1:nDFT)));
    elseif (i>nt-ceil(nDFT/2)+1)
        corr	= 1/(winCorr_fac*sum(window(1:nt+ceil(nDFT/2)-i)));
    else
        corr 	= 1;
    end
    for l=1:nModes
        a(:,l,i) = corr*winCorr_fac*dot(squeeze(P(:,l,:)),weight.*X_blk,2);
    end
    disp(['time ' num2str(i) '/' num2str(nt)])
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window] = hammwin(N)
%HAMMWIN Standard Hamming window of lenght N
    window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
end
