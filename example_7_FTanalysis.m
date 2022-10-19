%% EXAMPLE 7: Frequency-time analysis.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency–time analysis, low-rank
%         reconstruction and denoising of turbulent flows using SPOD,
%         Journal of Fluid Mechanics 926, A26, 2021
%
% A. Nekkanti (aknekkan@eng.ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 14-Oct-2022 (OTS)

clc, clear variables
addpath('utils')
disp('Loading the entire test database might take a second...')
load(fullfile('jet_data','jetLES.mat'),'p','p_mean','x','r','dt');

%% SPOD
%   We will use standard parameters: a Hamming window of length 256 and 50%
%   overlap.
nDFT  	= 128;
nOvlp 	= floor(nDFT/2);
weight 	= trapzWeightsPolar(r(:,1),x(1,:));

[L,P,f,~,A] ...
        = spod(p,nDFT,weight,nOvlp,dt);

figure 
loglog(f,L)
xlabel('frequency'), ylabel('SPOD mode energy')

%% Frequency-time analysis.
%   We will compute the mode expansion coefficients of the first 10 modes
%   using a windowing function and weights consistent with the SPOD.
nt              = size(p,1);
nModes          = 3;
a               = tcoeffs(p(1:nt,:,:),P,nDFT,weight,nModes);

%% Visualize the results.
t   = (1:nt)*dt;

figure

subplot(2,1,1)
pcolor(t,f,abs(squeeze(a(:,1,:)))); shading interp
daspect([100 1 1])
title(['frequency-time diagram (first mode, ' num2str(sum(L(:,1))/sum(L(:))*100,'%3.1f') '\% of energy)'])
xlabel('time'), ylabel('frequency'), caxis([0 0.75].*caxis)

subplot(2,1,2)
pcolor(t,f,squeeze(abs(sum(a,2)))); shading interp
daspect([100 1 1])
title(['frequency-time diagram (sum of first ' num2str(nModes) ' modes, ' num2str(sum(sum(L(:,1:nModes)))/sum(L(:))*100,'%3.1f') '\% of energy)'])
xlabel('time'), ylabel('frequency'), caxis([0 0.75].*caxis)

% %% Reconstruction from continuously-discrete temporal expansion coefficients
% %   The inverse SPOD using the continuously-discrete temporal expansion
% %   coefficients yields a rank-reduced data reconstruction. We use
% %   invspod() to reconstruct an entire block at a time, but only retain the
% %   central time instant. This is for demonstration purposes only and
% %   absolutely not recommended. Please refer to example 8 for the 'proper'
% %   way of reconstruction/filtering using the block-wise expansion
% %   coefficients.
% 
% nt      = 300;                          
% p_rec   = zeros([nt size(x)]);
% for t_i = 1:nt
%     p_block         = invspod(P,squeeze(a(:,:,t_i)),nDFT,nOvlp);
%     p_rec(t_i,:,:)  = p_block(ceil(nDFT/2)+1,:,:);
% end
% % 
% figure('name','Reconstruction from continuously-discrete temporal expansion coefficients')
% for t_i=1:nt
%     subplot(3,1,1)
%     pcolor(x,r,squeeze(p(t_i,:,:))-p_mean); shading interp, axis equal tight
%     if t_i==1; pmax = max(abs(caxis)); end
%     caxis(0.5*pmax*[-1 1]), colorbar  
% 
%     title('Original data')
%     caxis(0.5*pmax*[-1 1]), colorbar    
%     subplot(3,1,2)
%     pcolor(x,r,squeeze(p_rec(t_i,:,:))); shading interp, axis equal tight
%     title('Reconstructed data')
% 
%     caxis(0.5*pmax*[-1 1]), colorbar
%     subplot(3,1,3)
%     pcolor(x,r,(squeeze(p(t_i,:,:))-p_mean)-squeeze(p_rec(t_i,:,:))); shading interp, axis equal tight
%     title('Filtered/removed component')
%     caxis(0.5*pmax*[-1 1]), colorbar
%     drawnow
% end
