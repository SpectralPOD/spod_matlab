function [X] = invspod(P,A,window,nOvlp,dt)
%INVSPOD Inversion of SPOD using block-wise expansion coefficients
%
%   [X] = INVSPOD(P,A,WINDOW,NOVLP,DT) inverts the SPOD and returns the
%   original data X. P are the SPOD modes and A the block-wise SPOD
%   expansion coefficients returned by [L,P,F,Lc,A] = SPOD(X,...). WINDOW
%   is the data window and NOVLP the number of overlapping snapshots
%   between blocks. DT is the time step used in SPOD and defaults to 1.
%   The window-weighted average, equation (A1) in [1], is used in
%   overlapping regions. Multitaper estimation is not supported by this
%   version.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency–time analysis, low-rank 
%         reconstruction and denoising of turbulent flows using SPOD, 
%         Journal of Fluid Mechanics 926, A26, 2021
%
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 20-Aug-2026
%
% Revision history:
%   20-Aug-2026: Corrected overlap weighting when more than two blocks
%                contribute to a snapshot. We thank Joel Weightman for
%                finding and documenting this bug.
%   20-Aug-2026: Updated the inverse window and time-step scaling to match
%                the SPOD normalization introduced on 21-Aug-2025.

if nargin<5
    dt = 1;
end

dim     = size(P);
if ndims(A)==2
    nBlks   = 1;
else
    nBlks   = dim(end);
end
nFreqs  = dim(1);
nModes  = size(A,2);
if length(window)==1
    nDFT        = window;
    window      = hammwin(window);
else
    nDFT    = numel(window);
end
window  = window/sqrt(sum(window.^2));
nx      = prod(dim(2:end-1));
nt      = nDFT*nBlks-(nBlks-1)*nOvlp;
issymm  = nFreqs~=nDFT;

X       = zeros(nt,nx);
wGlob   = zeros(nt,1);
xHatBlk = zeros([nDFT nx]);
P       = permute(reshape(P,[nFreqs nx dim(end)]),[1 3 2]);

% loop over number of blocks and generate Fourier realizations
disp(' ')
disp('Reconstructing data from SPOD')
disp('------------------------------------')
for blk_i    = 1:nBlks
    % get time index for present block
    offset   = min((blk_i-1)*(nDFT-nOvlp)+nDFT,nt)-nDFT;
    timeIdx  = (1:nDFT) + offset;
    disp(['block ' num2str(blk_i) '/' num2str(nBlks) ' (' ...
        num2str(timeIdx(1)) ':' num2str(timeIdx(end)) ')'])
    xHatBlk(:)  = 0;
    
    % reconstruct Fourier realization
    xHatBlk(1:nFreqs,:)     = squeeze(sum(P(:,1:nModes,:).*squeeze(A(:,:,blk_i)),2));
    if issymm
        xHatBlk(nDFT/2+2:end,:)  = conj(xHatBlk(nDFT/2:-1:2,:));
    end
    
    % correction for windowing
    xBlk    = ifft(xHatBlk,nDFT,1)/sqrt(dt)./window;
    
    % weighted reconstruction in overlapping segments
    for ti_loc = 1:nDFT                             % local time index of current block
        ti_glob = timeIdx(ti_loc);                  % global time index
        w_curr  = window(ti_loc);                   % weight of current block
        if wGlob(ti_glob)~=0
            w_total     = wGlob(ti_glob)+w_curr;
            alpha       = w_curr/w_total;
            X(ti_glob,:)= alpha*xBlk(ti_loc,:) + (1-alpha)*X(ti_glob,:);
            wGlob(ti_glob) = w_total;
        else
            X(ti_glob,:)= xBlk(ti_loc,:);
            wGlob(ti_glob) = w_curr;
        end
    end
end

X   = reshape(X,[nt,dim(2:end-1)]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window] = hammwin(N)
%HAMMWIN Standard Hamming window of lenght N
    window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
end
