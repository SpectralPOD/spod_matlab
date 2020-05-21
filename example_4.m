%% EXAMPLE 4: Calculate the SPOD of large data and save results on hard drive.
%  The large-eddy simulation data provided along with this example is a
%  subset of the database of a Mach 0.9 turbulent jet described in [1] and 
%  was calculated using the unstructured flow solver Charles developed at 
%  Cascade Technologies. If you are using the database in your research or 
%  teaching, please include explicit mention of Brès et al. [1]. The test 
%  database consists of 5000 snapshots of the symmetric component (m=0) of 
%  a round turbulent jet. A physical interpretaion of the SPOD results is 
%  given in [2], and a comprehensive discussion and derivation of SPOD and
%  many of its properties can be found in [3].
%
%   References:
%     [1] G. A. Brès, P. Jordan, M. Le Rallic, V. Jaunet, A. V. G. 
%         Cavalieri, A. Towne, S. K. Lele, T. Colonius, O. T. Schmidt, 
%         Importance of the nozzle-exit boundary-layer state in subsonic 
%         turbulent jets, J. of Fluid Mech. 851, 83-124, 2018
%     [2] Schmidt, O. T. and Towne, A. and Rigas, G. and Colonius, T. and 
%         Bres, G. A., Spectral analysis of jet turbulence, J. of Fluid Mech. 855, 953–982, 2018
%     [3] Towne, A. and Schmidt, O. T. and Colonius, T., Spectral proper 
%         orthogonal decomposition and its relationship to dynamic mode
%         decomposition and resolvent analysis, J. of Fluid Mech. 847, 821–867, 2018
%
% O. T. Schmidt (oschmidt@ucsd.edu), A. Towne, T. Colonius
% Last revision: 20-May-2020

clc, clear variables
addpath('utils')
%   Note that we don't load the data 'p' itself.
load(fullfile('jet_data','jetLES.mat'),'p_mean','x','r','dt');

%% Memory-efficient SPOD version that stores and reloads FFT blocks from hard drive.
%   In this example, we use a function handle to provide individual
%   snapshots to SPOD(_). Additionally, we ask SPOD(_) to save the FFT
%   blocks on hard drive instead of keeping all data in memory. These
%   features enable the SPOD of very large data sets but require more
%   computing time and additional hard drive space. We reduce the
%   additional storage requirenment by saving only a few modes at selected
%   frequencies.

opts.savefft    = true;             % save FFT blocks insteasd of keeping them in memory
opts.deletefft  = false;            % keep FFT blocks in this example for demonstration purposes
opts.savedir    = 'results';        % save results to 'results' folder in the current directory
opts.savefreqs  = 10:5:20;          % save modes frequencies of indices [10 15 20]
opts.nt         = 2000;             % use 2000 snapshots using XFUN             
opts.mean       = p_mean;           % provide a long-time mean
opts.nsave      = 5;                % save the 5 most energetic modes

%   trapezoidal quadrature weights for cylindrical coordinates
intWeights = trapzWeightsPolar(r(:,1),x(1,:));

%   Use function handle to GETJET(_) to provide snapshots to SPOD(_). The
%   function file getjet.m can be found in the 'utils' folder. You can use 
%   getjet.m as a template to interface your own data with SPOD(_). A 
%   default (Hamming) window of length 256 with 128 snaphots overlap is
%   used in this example.
[L,P,f] = spod(@getjet,256,intWeights,128,dt,opts);

%% Plot the SPD spectrum and some modes as before.
%   Note that P is a function handle that loads the corresponding modes from
%   hard drive since we are in FFT saving mode (OPTS.savefft is true), and
%   that the spectrum is restricted to the 5 frequencies specified through
%   OPTS.savefreqs. 
figure 
loglog(f,L,'*')
xlabel('frequency'), ylabel('SPOD mode energy')

figure('name','Plot using PFUN function')
count = 1;
for fi = [10 15 20]
    for mi = [1 2]
        subplot(3,2,count)
        contourf(x,r,real(squeeze(P(fi,mi))),11,'edgecolor','none'), axis equal tight, caxis(max(abs(caxis))*[-1 1])
        xlabel('x'), ylabel('r'), title(['f=' num2str(f(fi),'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(L(fi,mi),'%.2g')])
        xlim([0 10]); ylim([0 2])
        count = count + 1;
    end
end

%% Repeat the same plot, but manually load the saved result files.
figure('name','Load SPOD modes from files')
count = 1;
for fi = [10 15 20]
    file = matfile(['results/nfft256_novlp128_nblks14/spod_f' num2str(fi,'%.4i')]);
    for mi = [1 2]
        subplot(3,2,count)
        contourf(x,r,real(squeeze(file.Psi(:,:,mi))),11,'edgecolor','none'), axis equal tight, caxis(max(abs(caxis))*[-1 1])
        xlabel('x'), ylabel('r'), title(['f=' num2str(f(fi),'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(L(fi,mi),'%.2g')])
        xlim([0 10]); ylim([0 2])
        count = count + 1;
    end
end

