function [DBL, SZ, NUP, T_DBL, T_SZ, T_NUP, T_NDN, MPn, MnP, MPt, MtP, P] = get_data(sim)
% GET_DATAV Loads data from cmpl_data and generates matrices for
% figures. It algo gives the correspondant meshgrids for surfing.
%   From 'sim' determines the data interval. Then determines the
%   parameter and fills the matrices with data from a loop.
%   Vertical interval: var_indep = 1
%   Horizontal interval: var_indep = -1
%   P term means parameter. Can take values of V or U.
clearvars -except sim
clc
%% Determines data interval
[P, F, L, chi_m, sentinel] = unify_data(sim);
t = 101; % Length of time array (not generalized)
n = 1:L; % Array of sites

%% Vertical interval (V)
if sentinel > 0
    l = length(P); % Length of parameter array
    % Initializacion of output matrices
    DBL = zeros(l,L);
    SZ = zeros(l,L);
    NUP = zeros(l,L);
    T_DBL = zeros(l,t);
    T_SZ = zeros(l,t);
    T_NUP = zeros(l,t);
    T_NDN = zeros(l,t);
    
    for i=1:l
        %sprintf('%d',i)
        fname = ['./output/data/data_L' num2str(L) '_J1_U'...
            num2str(F) '_V' num2str(P(i), '%.2f') '_chi' num2str(chi_m)...
            '.mat'];
        load (fname)
        DBL(i,:) = dbl;
        SZ(i,:) = sz;
        NUP(i,:) = nup;
        T_DBL(i,:) = real(TC_dbl);
        T_SZ(i,:) = real(TC_sz);
        T_NUP(i,:) = real(TC_nup);
        T_NDN(i,:) = real(TC_ndn);
    end
    % Meshgrids
    [MnP,MPn] = meshgrid(n,V);
    [MtP,MPt] = meshgrid(time,V);
    
%% Horizontal interval (U)
elseif sentinel < 0
    l = length(U);
    % Initializacion of output matrices
    DBL = zeros(l,L);
    SZ = zeros(l,L);
    NUP = zeros(l,L);
    T_DBL = zeros(l,t);
    T_SZ = zeros(l,t);
    T_NUP = zeros(l,t);
    T_NDN = zeros(l,t);
    
    for i=1:l
        fname = ['./output/cmpl_data/data_L' num2str(L) '_J1_U'...
            num2str(U(i), '%.2f') '_V' num2str(V) '_chi' num2str(chi_m)...
            '.mat'];
        load (fname)
        DBL(i,:) = dbl;
        SZ(i,:) = sz;
        NUP(i,:) = nup;
        T_DBL(i,:) = real(TC_dbl);
        T_SZ(i,:) = real(TC_sz);
        T_NUP(i,:) = real(TC_nup);
        T_NDN(i,:) = real(TC_ndn);
    end
    % Meshgrids
    [MnP,MPn] = meshgrid(n,U);
    [MtP,MPt] = meshgrid(time,U);
end
end
