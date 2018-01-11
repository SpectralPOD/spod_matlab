function [weight_rz] = trapzWeightsPolar(r,z)
%TRAPZWEIGHTSPOLAR Integration weight matrix for cylindical coordinates using trapazoidal rule
% OTS, 2015

nothetar = length(r);
weight_thetar = zeros(nothetar,1);
weight_thetar(1) = pi*( r(1) + (r(2)-r(1))/2)^2;
for i=2:nothetar-1
    weight_thetar(i) = pi*( r(i) + (r(i+1)-r(i))/2 )^2 - pi*( r(i) - (r(i)-r(i-1))/2 )^2;
end
weight_thetar(nothetar) = pi*r(end)^2 - pi*( r(end) - (r(end)-r(end-1))/2 )^2;

% dz
noz = length(z);
weight_z = zeros(noz,1);
weight_z(1) = (z(2)-z(1))/2;
for i=2:noz-1
    weight_z(i) = (z(i)-z(i-1))/2 + (z(i+1)-z(i))/2;
end
weight_z(noz) = (z(noz)-z(noz-1))/2;

weight_rz    = weight_thetar*weight_z';

