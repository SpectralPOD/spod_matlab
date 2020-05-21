%% EXAMPLE 1: Inspect data and plot SPOD spectrum.
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

% Clean up the worspace.
clc, clear variables

%% Load the test database.
%   'p' is the data matrix, 'x' and 'r' are the axial and radial 
%   coordinates, respectively.
load(fullfile('jet_data','jetLES.mat'),'p','x','r');

%% Inspect the database.
%   The jet is resolvend by 39 points in the radial and by 175 points in 
%   the axial direction. To use the SPOD(_) function, it is important that
%   the first dimension of the data is time, i.e. 5000 snapshots in this 
%   example. If your data is sorted differently, please use PERMUTE() to 
%   swap indices.
size(p)

%% Animate the first 100 snapshots of the pressure field.
figure('name','Pressure of the symmetric component of a turbulent jet')
for ti=1:100
    pcolor(x,r,squeeze(p(ti,:,:)))
    axis equal tight, shading interp, caxis([4.43 4.48])
    xlabel('x'), ylabel('r')
    pause(0.05)
    drawnow
end

%% Calculate the SPOD spectrum of the data.
%   Check the output of SPOD(_) in the Command Window after execution: the
%   routine has segmented the data into 38 blocks of 256 snapshots each.
%   The segments, or blocks, overlap bei 50%, i.e. 128 snaphots overlap. 
[L]     = spod(p);

figure 
loglog(L)
xlabel('frequency index'), ylabel('SPOD mode energy')
