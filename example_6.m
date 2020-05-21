%% EXAMPLE 6: Calculate and plot confidence intervals for SPOD eigenvalues.
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
disp('Loading the entire test database might take a second...')
load(fullfile('jet_data','jetLES.mat'),'p','x','r','dt');

%% SPOD with 99% confidence interval.
%   This example is similar to example 3, but we also obtain the 99%
%   confidence interval of the SPOD eigenvalues by calling
%   [L,P,f,Lc] = SPOD(_) and specifying OPTS.conflvl.

opts.conflvl    = 0.99;             % return 99% confidence interval

%   trapezoidal quadrature weights for cylindrical coordinates
intWeights = trapzWeightsPolar(r(:,1),x(1,:));

%   SPOD using a retangular window of length 256 and 50 snaphots overlap
[L,P,f,Lc] = spod(p,ones(1,256),intWeights,50,dt,opts);

%% Plot the SPOD spectrum with upper and lower 99% confidence levels for every 5th mode.
figure
for mi = 1:5:size(L,2)
    lh = loglog(f,L(:,mi),'LineWidth',1); hold on
    loglog(f,Lc(:,mi,1),'LineWidth',0.1,'Color',get(lh,'Color'),'LineStyle','--'); % lower confidence level
    loglog(f,Lc(:,mi,2),'LineWidth',0.1,'Color',get(lh,'Color'),'LineStyle','--'); % upper confidence level
end
set(gca,'XScale','log','YScale','log');
xlabel('frequency'), ylabel('SPOD mode energy')
