%% EXAMPLE 3: Manually specify spectral estimation parameters and use cell volume weighted inner product.
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
%         turbulent jets, submitted to JFM, 2017
%     [2] Schmidt, O. T. and Towne, A. and Rigas, G. and Colonius, T. and 
%         Bres, G. A., Spectral analysis of jet turbulence, 
%         arXiv:1711.06296, 2017
%     [3] Towne, A. and Schmidt, O. T. and Colonius, T., Spectral proper 
%         orthogonal decomposition and its relationship to dynamic mode
%         decomposition and resolvent analysis, arXiv:1708.04393, 2017
%
% O. T. Schmidt (oschmidt@caltech.edu), A. Towne, T. Colonius
% Last revision: 5-Jan-2018

clc, clear variables
addpath('utils')
disp('Loading the entire test database might take a second...')
load(fullfile('jet_data','jetLES.mat'),'p','x','r','dt');

%% SPOD of the test database.
%   In this example, we manually specify a rectangular window of length 256 
%   and an overlap of 50 snaphots. Furthermore, we use trapezoidal
%   quadrature weights to define a physical inner product corresponding to
%   the volume integral over the pressure squared.

%   trapezoidal quadrature weights for cylindrical coordinates
intWeights = trapzWeightsPolar(r(:,1),x(1,:));

%   SPOD using a retangular window of length 256 and 50 snaphots overlap
[L,P,f] = spod(p,ones(1,256),intWeights,50,dt);

%% Plot the SPOD spectrum and some modes as before.
figure 
loglog(f,L)
xlabel('frequency'), ylabel('SPOD mode energy')

figure
count = 1;
for fi = [10 15 25]
    for mi = [1 2]
        subplot(3,2,count)
        get(gca,'Position')
        contourf(x,r,real(squeeze(P(fi,:,:,mi))),11,'edgecolor','none'), axis equal tight, caxis(max(abs(caxis))*[-1 1])
        xlabel('x'), ylabel('r'), title(['f=' num2str(f(fi),'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(L(fi,mi),'%.2g')])
        xlim([0 10]); ylim([0 2])
        count = count + 1;
    end
end