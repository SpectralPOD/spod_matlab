%% EXAMPLE 9: Multitaper-Welch.
%  Multitaper-Welch estimators provide lower variance estimates at a fixed
%  frequency resolution or higher frequency resolution at similar variance
%  compared to the standard algorithm. In this example, we retain the high
%  frequency resolution of a three block Welch estimate but significantly
%  reduce the variance of the SPOD spectrum by using 10 Slepian tapers.
%
%   References:
%     [1] O. T. Schmidt, Spectral proper orthogonal decomposition using
%         multitaper estimates, Theor. Comput. Fluid Dyn., 2022, 1-14, 
%         DOI 10.1007/s00162-022-00626-x, https://rdcu.be/cUtP3
%
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 5-Sep-2022

clc, clear variables
addpath('utils')
disp('Loading the entire test database might take a second...')
load(fullfile('jet_data','jetLES.mat'),'p','x','r','dt');

%   trapezoidal quadrature weights for cylindrical coordinates
intWeights = trapzWeightsPolar(r(:,1),x(1,:));

%% Standard SPOD
%   SPOD with a large block size to get a high frequency resolution and
%   resolve the low-frequency regime.
nFFT    = 2048;
nOvlp   = nFFT/2;
[L,P,f] = spod(p,nFFT,intWeights,nOvlp,dt);

%   Plot the SPOD spectrum and three leading modes at the frequency of
%   interest.
f_plot  = 0.24;
[~,fi]  = min(abs(f-f_plot));
nBlk    = size(L,2);

figure
subplot(5,2,[1 3])
loglog(f,L)
xlim([f(2) f(end)]), ylim([1e-9 1e-3])
xlabel('frequency'), ylabel('SPOD mode energy')
title('Welch (standard)')
hold on
plot([f(fi) f(fi)],ylim,'k:')

count   = 5;
for mi  = 1:3
    subplot(5,2,count)
    contourf(x,r,real(squeeze(P(fi,:,:,mi))),11,'edgecolor','none'), axis equal tight, caxis(max(abs(caxis))*[-1 1])
    xlabel('x'), ylabel('r'), title(['f=' num2str(f(fi),'%.2f') ', mode ' num2str(mi)])
    xlim([0 10]); ylim([0 2])
    count = count + 2;
end
drawnow

%% Multitaper-Welch SPOD
%   SPOD using a retangular window of length 256 and 50 snaphots overlap
%   10 Slapian tapers by setting the time-halfbandwidth product to 5.5.
bw      = 5.5;
[L,P,f] = spod(p,[nFFT bw],intWeights,nOvlp,dt);

%   Plot the SPOD spectrum and modes as before. Compared to the standard
%   algorithm, the variance of the spectrum has been reduced significanlty
%   and the modes are better converged.
subplot(5,2,[2 4])
loglog(f,L(:,1:nBlk)), hold on
loglog(f,L(:,nBlk+1:end),'Color',[0.75 0.75 0.75]), hold on
xlim([f(2) f(end)]), ylim([1e-9 1e-3])
xlabel('frequency'), ylabel('SPOD mode energy')
title(['Multitaper-Welch, b_w=' num2str(bw)])

hold on
plot([f(fi) f(fi)],ylim,'k:')

count = 6;
for mi = 1:3
    subplot(5,2,count)
    contourf(x,r,real(squeeze(P(fi,:,:,mi))),11,'edgecolor','none'), axis equal tight, caxis(max(abs(caxis))*[-1 1])
    xlabel('x'), ylabel('r'), title(['f=' num2str(f(fi),'%.2f') ', mode ' num2str(mi)])
    xlim([0 10]); ylim([0 2])
    count = count + 2;
end
