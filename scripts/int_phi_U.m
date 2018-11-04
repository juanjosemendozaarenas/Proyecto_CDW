clear
clc
close all
V = 4;
chi_m = 1000;
U = 7.5:0.05:8.5;
l = length(U);
S_dbl = zeros(l,101); %Generalizar!
S_sz = zeros(l,101); %Generalizar!

for i = 1:length(U)
    fname = ['output/cmpl_data/data_L64_J1_U'...
        num2str(U(i), '%.2f') '_V' num2str(V)...
        '_chi' num2str(chi_m) '.mat'];
    load(fname,'TC_dbl','TC_sz')
    if i==1
        load(fname,'time')
    end
    phi_dbl = sum_rec(real(TC_dbl));
    phi_sz = sum_rec(real(TC_sz));
    S_dbl(i,:) = phi_dbl;
    S_sz(i,:) = phi_sz;
end
[mU,mT] = meshgrid(time,U);
figure(1)
surf(mT,mU,S_dbl)
title('\phi para doble ocupación V=4')
xlabel('U')
ylabel('Tiempo')

figure(2)
surf(mT,mU,S_sz)
title('\phi para magnetización V=4')
xlabel('U')
ylabel('Tiempo')