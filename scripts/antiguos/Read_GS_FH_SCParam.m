function [O_s,O_snn,O_tnn] = Read_GS_FH_SCParam(V)
% d and m are order parameters (double ocup & Magnetization)
% for a given U. E0 is the ground state energy.
% These 3 values are intialized as NaNs to prevent errors if the
% initfile is not found
O_s=NaN;
O_snn=NaN;
O_tnn=NaN;
%% File parameters
L = 32;
J = 1;
% V = 6;
U = 1;
chi_max = 1500;

fname = ['/output/L' num2str(L) '_J' num2str(J) ...
    '_chi' num2str(chi_max) '_4-mar-18/GS_FH_L' num2str(L) '_J' num2str(J)...
    '_U' num2str(U) '_V' num2str(V) '_chi' num2str(chi_max) '.mat'];
check = exist(fname);
if check == 2
    load(fname);
    %% O_all sites
    O_s=sum(Osos(:))/L;
    
    %% O_nn singlet
    O_snn=(sum(nupndn1)-sum(Onn2)-sum(Onn3)+sum(ndnnup4))/L;

     %% O_nn triplet
    O_tnn=(sum(nupndn1)+sum(Onn2)+sum(Onn3)+sum(ndnnup4))/L;
else
    sprintf('File not found')
end

