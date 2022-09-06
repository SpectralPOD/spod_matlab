%% EXAMPLE 8: Band-pass filtering.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency–time analysis, low-rank reconstruction and denoising of turbulent flows using SPOD, 
%         Journal of Fluid Mechanics 926, A26, 2021
%
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 5-Sep-2022

clc, clear variables
addpath('utils')
disp('Loading the entire test database might take a second...')
load(fullfile('jet_data','jetLES.mat'),'p','p_mean','x','r','dt');

%   trapezoidal quadrature weights for cylindrical coordinates
intWeights      = trapzWeightsPolar(r(:,1),x(1,:));

%% SPOD
%   We will use standard parameters: a Hann window of length 256 and 50%
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
p_rec           = invspod(P,A,nDFT,nOvlp);

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

