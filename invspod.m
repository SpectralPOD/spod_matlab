function [X] = invspod(P,A,window,nOvlp)
%INVSPOD Inversion of SPOD using block-wise expansion coefficients
%
%   [X] = INVSPOD(P,A,WINDOW,NOVLP) inverts the SPOD and returns the
%   original data X. P are the SPOD modes and A the block-wise SPOD
%   expansion coefficients returned by [L,P,F,Lc,A] = SPOD(X,...). WINDOW
%   is the data window and NOVLP the number of overlapping snapshots
%   between blocks. The window-weighted average, equation (A1) in [1], is
%   used in overlapping regions. Multitaper estimation is not supported by
%   this version.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency–time analysis, low-rank 
%         reconstruction and denoising of turbulent flows using SPOD, 
%         Journal of Fluid Mechanics 926, A26, 2021
%
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 7-Oct-2022 

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
nx      = prod(dim(2:end-1));
nt      = nDFT*nBlks-(nBlks-1)*nOvlp;
issymm  = nFreqs~=nDFT;

X       = zeros(nt,nx);
xHatBlk = zeros([nDFT nx]);
P       = permute(reshape(P,[nFreqs nx dim(end)]),[1 3 2]);
timeIdx_prev    ...
        = [];

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
    xHatBlk = xHatBlk*sum(window);
    xBlk    = ifft(xHatBlk,nDFT,1)./window;
    
    % weighted reconstruction in overlapping segments
    for ti_loc = 1:nDFT                             % local time index of current block
        ti_glob = timeIdx(ti_loc);                  % global time index
        ti_prev = find(timeIdx_prev==ti_glob,1);    % local time index of previous block
        if ~isempty(ti_prev)
            w_prev      = window(ti_prev);
            w_curr      = window(ti_loc);
            alpha       = w_curr/(w_curr+w_prev);
            X(ti_glob,:)= alpha*xBlk(ti_loc,:) + (1-alpha)*X(ti_glob,:);
        else
            X(ti_glob,:)= xBlk(ti_loc,:);
        end
    end
    timeIdx_prev = timeIdx;
end

X   = reshape(X,[nt,dim(2:end-1)]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window] = hammwin(N)
%HAMMWIN Standard Hamming window of lenght N
    window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
end