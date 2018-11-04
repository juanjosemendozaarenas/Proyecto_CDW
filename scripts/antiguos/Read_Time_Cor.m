function [tc] = Read_Time_Cor( V, U, t)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
L = 64;
J = 1;
chi = 1500; % Máximo chi para DMRG
chi_max = 1200; % Máximo chi para evolución temporal

% Archivo corte vertical
fname = ['output/TE1a/TimeCorrel_FH1a_z_L' num2str(L) '_J'...
    num2str(J) '_U' num2str(U) '_V' num2str(V,'%0.2f') '_chi' ...
    num2str(chi_max) '_dt0.01_era0.mat'];

% Archivo corte horizontal
% fname = ['output/TE3a/TimeCorrel_FH3a_z_L' num2str(L) '_J'...
%     num2str(J) '_U' num2str(U,'%0.2f') '_V' num2str(V) ...
%     '_chi' num2str(chi_max) '_dt0.01_era0.mat'];
check = exist(fname);
if check == 2
    load(fname);
    tc = real(TimeCorrel);
    
else
    tc = zeros(1,t);
    warning = ['File for U=' num2str(U) ' and V=' num2str(V)...
        ' not found. Taking zero por time correlation value'];
    sprintf(warning)

end

