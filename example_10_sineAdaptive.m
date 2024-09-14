%% EXAMPLE 10: Adaptive sine-taper SPOD.
%  The time-resolved particle image velocimetry data provided along with
%  this example are a subset of the database of a Mach 0.6 turbulent open
%  cavity flow by Zhang et al. [1]. If you are using the database in your
%  research or teaching, please include explicit mention of Zhang et al.
%  [1]. The test database consists of 4000 snapshots of the streamwise
%  velocity component of a turbulent open cavity flow. This is a canonical
%  example of a mixed broadband-tonal turbulent flow, which requires sharp
%  peak resolution at tonal frequencies and smooth spectrum estimates of
%  broadband regions. If you are using this code for your research, please
%  cite Yeung & Schmidt [2].
%
%   References:
%     [1] Zhang, Y., Cattafesta, L. N., and Ukeiley, L.,  Spectral analysis
%         modal methods (SAMMs) using non-time-resolved PIV, Exp. Fluids 61, 226, 1–12, 2020,
%         DOI 10.1007/s00348-020-03057-8
%     [2] Yeung, B. C. Y., Schmidt, O. T. Adaptive spectral proper orthogonal decomposition
%         of broadband-tonal flows, Theor. Comput. Fluid Dyn. 38, 355–374, 2024,
%         DOI 10.1007/s00162-024-00695-0
%
% B. Yeung (byeung@ucsd.edu) and O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 10-Sep-2024

clc, clear variables
addpath('utils')
load(fullfile('cavity_data','cavityPIV.mat'),'u','x','y','dt');

%% Standard Welch SPOD of the test database using one sine window.
%   We manually specify a window length of 256 and an overlap
%   of 128 snaphots. A sine function is automatically used as the window.

%   SPOD using a sine window of length 256 and 128 snaphots overlap
nDFT    = 256;
nOvlp   = 128;
[L,P,f] = spod_adapt(u,nDFT,[],nOvlp,dt);

%% Plot the SPOD spectrum and some modes.
figure 
loglog(f,L)
xlabel('frequency'), ylabel('SPOD mode energy'), title('Welch SPOD')

%% Multitaper SPOD using six sine windows.
%   We can achieve higher frequency resolution and resolve down to lower
%   minimum frequencies using six sine windows of length equal to the
%   number of snapshots. Specify nWin as a vector of length 4000.
%   Multitaper SPOD neither requires nor benefits from overlap, so we leave
%   nOvlp empty.

%   SPOD using six sine windows of length 4000 and 0 overlap
nWin    = 6*ones(size(u,1),1);
[L,P,f] = spod_adapt(u,nWin,[],[],dt);

%% Plot the SPOD spectrum.
figure 
loglog(f,L)
xlabel('frequency'), ylabel('SPOD mode energy'), title('multitaper SPOD')

%% Adaptive multitaper SPOD.
%   To combine both sharp peak resolution and smooth spectrum estimates of
%   broadband regions, we allow the window number to vary with frequency.
%   The adaptive algorithm requires only a tolerance and determines the
%   window number autonomously. No spectral estimation parameters
%   are needed.

%   SPOD using the adaptive algorithm with a tolerance of 1e-4
%   Can take a long time to compute
clear opts
opts.adaptive   = true;
opts.tol        = 1e-4;
[L,P,f,~,nWin]  = spod_adapt(u,[],[],[],dt,opts);

%% Plot the SPOD spectrum and window number.
figure
yyaxis left
loglog(f,L)
xlabel('frequency'), ylabel('SPOD mode energy'), title('adaptive multitaper SPOD')

yyaxis right
semilogx(f,nWin)
ylabel('number of windows')

%% Plot some SPOD modes (Rossiter modes)
figure
count = 1;
for fi = [470 770 1063]
    for mi = [1 2]
        subplot(3,2,count)
        contourf(x,y,real(squeeze(P(fi,:,:,mi))),11,'edgecolor','none'), axis equal tight, clim(max(abs(clim))*[-1 1])
        xlabel('x'), ylabel('y'), title(['f=' num2str(f(fi),'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(L(fi,mi),'%.2g')])
        count = count + 1;
    end
end

% %% Multitple SPOD using non-uniform window number.
% %   If the (potentially non-uniform) window number is known a priori, e.g.,
% %   based on physical knowledge, it can be provided directly. As an example,
% %   we use the nWin computed previously to recover the same results as before.
% 
% %   SPOD using sine windows of length 4000, 0 overlap, and
% %   frequency-dependent window number
% if length(nWin)<size(u,1), nWin(size(u,1)) = 0; end % zero-pad nWin to length 4000
% [L,P,f] = spod_adapt(u,nWin,[],[],dt);
% 
% %% Plot the SPOD spectrum.
% %   The spectrum should be identical to the previous figure.
% figure 
% loglog(f,L)
% xlabel('frequency'), ylabel('SPOD mode energy'), title('multitaper SPOD with non-uniform window number')
