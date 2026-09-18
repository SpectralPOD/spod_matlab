%% EXAMPLE 8: Band-pass filtering.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency–time analysis, low-rank reconstruction and denoising of turbulent flows using SPOD, 
%         Journal of Fluid Mechanics 926, A26, 2021
%
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 20-Aug-2026
%
% Revision history:
%   20-Aug-2026: Added a reconstruction test for overlaps above 50% after
%                a bug report by Joel Weightman.
%   20-Aug-2026: Updated inverse-SPOD calls for corrected window and DT
%                normalization.

clc, clear variables
addpath('utils')

%% Validate reconstruction for multiple overlaps
%   Full-rank inverse SPOD must reproduce the mean-subtracted input for any
%   valid overlap. Overlaps above 50% cause three or more blocks to
%   contribute to some snapshots and therefore exercise the corrected
%   overlap-weight accumulation in INVSPOD.
rng(1)
nDFT_test      = 64;
nOvlp_test     = [nDFT_test/2 3*nDFT_test/4 7*nDFT_test/8];
dt_test        = 0.2;
x_test         = randn(512,64);
recError       = zeros(size(nOvlp_test));
for iOvlp = 1:numel(nOvlp_test)
    [~,P_test,~,~,A_test] = spod(x_test,nDFT_test,[],nOvlp_test(iOvlp),dt_test);
    x_rec_test             = invspod(P_test,A_test,nDFT_test,nOvlp_test(iOvlp),dt_test);
    x_ref_test             = x_test(1:size(x_rec_test,1),:)-mean(x_test,1);
    recError(iOvlp)        = norm(x_rec_test-x_ref_test,'fro')/norm(x_ref_test,'fro');
end

fprintf('\nInverse-SPOD reconstruction check:\n')
fprintf('  overlap = %5.1f%%, relative error = %.3e\n', ...
    [100*nOvlp_test/nDFT_test; recError])
assert(all(recError<1e-10),'INVSPOD reconstruction test failed.')

%% Jet-data band-pass filtering
disp('Loading the entire test database might take a second...')
load(fullfile('jet_data','jetLES.mat'),'p','p_mean','x','r','dt');

%   trapezoidal quadrature weights for cylindrical coordinates
intWeights      = trapzWeightsPolar(r(:,1),x(1,:));

%% SPOD
%   We will use standard parameters: a Hamming window of length 256 and 50%
%   overlap.
nDFT            = 256;
nOvlp           = nDFT/2;
[L,P,f,~,A]     = spod(p,nDFT,intWeights,nOvlp,dt);

figure
loglog(f,L)
title('SPOD of full data')
xlabel('frequency'), ylabel('SPOD mode energy')
ylims           = ylim;

%% Filtering
%   An band-pass filter that focusses on the 'low-rank' portion of the
%   spectrum is implemented by setting to zero the SPOD expansion
%   coefficients in the range f<=0.08 and f>=1.0.
f_lowpass       = 1;
f_highpass      = 0.08;  
A(f>=f_lowpass|f<=f_highpass,:,:) ...
                = 0;

%   The inverse SPOD using the modified SPOD expansion coefficients yields
%   the band-pass filtered data.
nt              = size(p,1);
p_rec           = invspod(P,A,nDFT,nOvlp,dt);

%% Animate
%   Animate the original, filtered, and removed data.
figure
for t_i=1:1:30
    subplot(3,1,1)
    pcolor(x,r,squeeze(p(t_i,:,:))-p_mean); shading interp, axis equal tight
    title('Original data')
    if t_i==1; pmax = max(abs(caxis)); end, caxis(0.5*pmax*[-1 1]), colorbar    
    subplot(3,1,2)
    pcolor(x,r,squeeze(p_rec(t_i,:,:))); shading interp, axis equal tight
    title('Reconstructed data')
    caxis(0.5*pmax*[-1 1]), colorbar
    subplot(3,1,3)
    pcolor(x,r,(squeeze(p(t_i,:,:))-p_mean)-squeeze(p_rec(t_i,:,:))); shading interp, axis equal tight
    title('Filtered/removed component')
    caxis(0.5*pmax*[-1 1]), colorbar
    drawnow
end

%% SPOD of filtered data
%   The effect of the band-pass filter should to a large degree remove
%   high-frequency components above the filter cut-on frequency.
[L,P,f,~,A]     = spod(p_rec,nDFT,intWeights,nOvlp,dt);

figure
loglog([f_lowpass f_lowpass],ylims,'k--'); hold on
loglog([f_highpass f_highpass],ylims,'k:');
loglog(f,L)
title('SPOD of filtered data')
xlabel('frequency'), ylabel('SPOD mode energy')
ylim(ylims);
legend('f_{low-pass}','f_{high-pass}')
