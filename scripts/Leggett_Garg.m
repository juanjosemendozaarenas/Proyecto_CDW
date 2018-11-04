function [LG, MtP, MPt] = Leggett_Garg( TC, P)
% LEGGET_GARG calculates the Legget-Garg inequalities for a given plane
% of Time Correlations (TC) and a parameter(P). Two meshgrids are given
% as well.
%   TC is a plane of parameter P times time T. The inequalities
%   go from tau=0 to tau=T/2.

[f,c] = size(TC); %Dimension of TC plane (parameter,time)
t = ceil(c/2); % The dimension of the LGI plane
tau = 0:0.02:1; % Tau array (only for meshgrid)
[MtP,MPt] = meshgrid(tau,P); % Meshgrid
LG = zeros(f,t); % Initialization of LGI plane

for i=1:f % loop over parameter
    for j=1:t % loop over time(tau)
        LG(i,j) = TC(i,(2*j-1))-2*TC(i,j);
    end
end

end

