clear all;
close all;
clc;
%%
% %1a
% V = -1:0.02:-0.5;U = -2;

% %1b
% V = -0.1:0.01:0.1;U = -2;
%ext
V = [-0.5:0.05:-0.15 -0.1:0.01:0.1 0.15:0.05:0.5]; U = 1;

% %2a
% V = -1.6:0.01:-1.4; U=1;
% %ext
% V = -1.6:0.01:-1.31; U=1;

% %2b
% V = -1.3:0.02:-0.7;U = 1;
% %ext
% V = -1.3:0.02:-0.4;U = 1;

% %2d
% V = -0.37:0.03:0.17; U = 1;

% %2c
% V = 0.2:0.02:0.7;
% %ext
% V = 0.2:0.02:0.9;

v = length(V);
L = 64;
N = 1:1:L;
DBL = zeros(v,L);
NUP = zeros(v,L);
SZ = zeros(v,L);
[mN,mV] = meshgrid(N,V);


for j=1:v
    fname = ['output/ExOps2b/GS_FH2b_L64_J1_U'...
        num2str(U) '_V' num2str(V(j),'%0.2f') '_chi1500.mat'];
    check = exist(fname);
    if check == 2
        load(fname)
        DBL(j,:) = dbl;
        NUP(j,:) = nup;
        SZ(j,:) = sz;
    else
        warning = ['File for U=' num2str(U) ' and V='...
            num2str(V(j),'%0.2f') ' not found.'];
    sprintf(warning)
    end
end

%% Figure 1
chosen_n = [30,25,64];
figure(1)
surf(mV,mN,DBL)
title('Doble ocupación')
xlabel('V')
ylabel('N')

figure(2)
plot(V,DBL(:,chosen_n(1)),'linewidth', 2);
title(['Doble ocupación en N = ' num2str(chosen_n(1))])
xlabel('V')
ylabel('<dbl>')

%% Figure 2
figure(3)
surf(mV,mN,NUP)
title('Espines arriba')
xlabel('V')
ylabel('N')

figure(4)
plot(V,NUP(:,chosen_n(2)),'linewidth', 2);
title(['Espines arriba en N = ' num2str(chosen_n(2))])
xlabel('V')
ylabel('<nup>')

%% Figure 3
figure(5)
surf(mV,mN,SZ)
title('Magnetización')
xlabel('V')
ylabel('N')

figure(6)
plot(V,SZ(:,chosen_n(3)),'linewidth', 2);
title(['Magnetización en N = ' num2str(chosen_n(3))])
xlabel('V')
ylabel('<sz>')