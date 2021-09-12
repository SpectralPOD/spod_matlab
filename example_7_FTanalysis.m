%% EXAMPLE 7: Frequency-time analysis.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency–time analysis, low-rank reconstruction and denoising of turbulent flows using SPOD, 
%         Journal of Fluid Mechanics 926, A26, 2021
%
% A. Nekkanti (aknekkan@eng.ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 12-Sep-2021

clc, clear variables
addpath('utils')
disp('Loading the entire test database might take a second...')
load(fullfile('jet_data','jetLES.mat'),'p','x','r','dt');

%% SPOD
%   We will use standard parameters: a Hann window of length 256 and 50%
%   overlap.
nFFT  	= 128;
nOvlp 	= floor(nFFT/2);
weight 	= trapzWeightsPolar(r(:,1),x(1,:));
window	= hann(nFFT);

[L,P,f] = spod(p,window,weight,nOvlp,dt);

figure 
loglog(f,L)
xlabel('frequency'), ylabel('SPOD mode energy')

%% Frequency-time analysis.
%   We will compute the mode expansion coefficients of the first 3 modes
%   using a windowing function and weights consistent with the SPOD.
nModes          = 3;
a               = tcoeffs(p,P,window,weight,nModes);

%% Visualize the results.
nt  = size(p,1);
t   = (1:nt)*dt;

figure

subplot(2,1,1)
pcolor(t,f,abs(squeeze(a(:,:,1)))); shading interp
daspect([100 1 1])
title(['frequency-time diagram (first mode, ' num2str(sum(L(:,1))/sum(L(:))*100,'%3.1f') '% of energy)'])
xlabel('time'), ylabel('frequency'), caxis([0 0.75].*caxis)

subplot(2,1,2)
pcolor(t,f,squeeze(abs(sum(a,3)))); shading interp
daspect([100 1 1])
title(['frequency-time diagram (sum of first ' num2str(nModes) ' modes, ' num2str(sum(sum(L(:,1:3)))/sum(L(:))*100,'%3.1f') '% of energy)'])
xlabel('time'), ylabel('frequency'), caxis([0 0.75].*caxis)

%% Alternatively, plot for pre-multiplied spectra
% pcolor(t,f,abs(squeeze(a(:,:,1)).*f')); shading interp
% set(gca,'YScale','log')
% daspect([50 1 1])
% title('premultiplied frequency-time diagram')
