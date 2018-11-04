clear
clc
close all
U = -2;
%U = 1;
chi_m = 1500;
V = [-1:0.02:-0.5 -0.45:0.05:-0.15 -0.1:0.01:0.1 0.15:0.05:0.5];
%V = [-2:0.02:-1.62 -1.6:0.01:-1.31 -1.3:0.02:-0.4 -0.37:0.03:0.17 0.2:0.02:0.9];
l = length(V);
S_dbl = zeros(l,101); %Generalizar!
S_sz = zeros(l,101); %Generalizar!

for i = 1:length(V)
    fname = ['output/cmpl_data/data_L64_J1_U' num2str(U)...
        '_V' num2str(V(i), '%.2f') '_chi' num2str(chi_m) '.mat'];
    load(fname,'TC_dbl','TC_sz')
    if i==1
        load(fname,'time')
    end
    phi_dbl = sum_rec(real(TC_dbl));
    phi_sz = sum_rec(real(TC_sz));
    S_dbl(i,:) = phi_dbl;
    S_sz(i,:) = phi_sz;
end
[mV,mT] = meshgrid(time,V);
figure(1)
surf(mT,mV,S_dbl)
title(['\phi para doble ocupación U=' num2str(U)])
xlabel('V')
ylabel('Tiempo')

figure(2)
surf(mT,mV,S_sz)
title(['\phi para magnetización U=' num2str(U)])
xlabel('V')
ylabel('Tiempo')

    