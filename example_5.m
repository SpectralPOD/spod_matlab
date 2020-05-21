%% EXAMPLE 5: Calculate full SPOD spectrum of large data.
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
load(fullfile('jet_data','jetLES.mat'),'p_mean','x','r','dt');

%% Memory-efficient SPOD version that stores and reloads FFT blocks from hard drive.
%   This example is almost identical to the previous one, but we're solely
%   interested in obtaining a full SPOD spectrum without storing any modes.
%   We therefore do not restrict the frequency range by specifying 
%   OPTS.savefreqs, and we set OPTS.save to null.

opts.savefft    = true;             % save FFT blocks insteasd of keeping them in memory
opts.savedir    = 'results';        % save results to 'results' folder in the current directory
opts.nt         = 2000;             % use 2000 snapshots using XFUN             
opts.mean       = p_mean;           % provide a long-time mean
opts.nsave      = 0;                % do not save the modes; we're only interested in the full spectrum

%   trapezoidal quadrature weights for cylindrical coordinates
intWeights = trapzWeightsPolar(r(:,1),x(1,:));

%   Calculate SPOD with default window of length 128 and with 50% overlap
[L,~,f] = spod(@getjet,128,intWeights,64,dt,opts);

%% Plot the SPD spectrum and some modes as before.
%   Note that P is a function handle that loads the corresponding modes from
%   hard drive since we are in FFT saving mode (OPTS.savefft is true).
figure 
loglog(f,L)
xlabel('frequency'), ylabel('SPOD mode energy')
