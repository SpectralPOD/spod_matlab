%% EXAMPLE 2: Plot SPOD spectrum and inspect SPOD modes.
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
% Last revision: 5-Sep-2022 (OTS)

clc, clear variables
addpath('utils')
load(fullfile('jet_data','jetLES.mat'),'p','x','r','dt');

%% SPOD of the test database.
%   Calculate the SPOD of the data matrix 'p' and use the timestep 'dt'
%   between snapshots to obtain the physical frequency 'f'. 'L' is the
%   matrix of modal energies, as before, and 'P' the data matrix of SPOD
%   modes. We leave all other options empty for now.
[L,P,f] = spod(p,[],[],[],dt);

%   First, we plot the SPOD spectrum.
figure
loglog(f,L)
xlabel('frequency'), ylabel('SPOD mode energy')

%   Second, we visualize the 1st and 2nd SPOD modes at three frequencies.
figure
count = 1;
for fi = [10 15 20]
    for mi = [1 2]
        subplot(3,2,count)
        contourf(x,r,real(squeeze(P(fi,:,:,mi))),11,'edgecolor','none'), axis equal tight, caxis(max(abs(caxis))*[-1 1])
        xlabel('x'), ylabel('r'), title(['f=' num2str(f(fi),'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(L(fi,mi),'%.2g')])
        xlim([0 10]); ylim([0 2])
        count = count + 1;
    end
end

%% Animate the same modes.
%   Note how all wavepackets travel at approximately the same phase
%   speed c_ph. The reason is that their streamwise wavenumber k_x changes 
%   with frequency such that c_ph = omega/k_x is approximately constant.
figure
nt      = 30;
T       = 1/f(10);              % period of the 10th frequency
time    = linspace(0,T,nt);     % animate over one period
count = 1;
for ti = 1:nt
    for fi = [10 15 20]
        for mi = [1 2]
            subplot(3,2,count)
            pcolor(x,r,real(squeeze(P(fi,:,:,mi)*exp(2i*pi*f(fi)*time(ti))))), shading interp, axis equal tight, caxis(max(abs(caxis))*[-1 1])
            xlabel('x'), ylabel('r'), title(['f=' num2str(f(fi),'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(L(fi,mi),'%.2g')])
            xlim([0 10]); ylim([0 2])
            count = count + 1;
            hold on
        end
    end
    drawnow
    hold off
    count = 1;
end
