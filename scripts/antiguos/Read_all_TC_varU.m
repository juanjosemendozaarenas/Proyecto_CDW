clear
close all
clc
%% System Parameters
V = 4;
U = 7.5:0.05:8.5;
deltaH = U(2)-U(1);
l = length(U);
deltaT = 0.02;
T = 0:deltaT:2;
t = length(T);
%% Reading files to set TC matrix
plano = [];
j=1;
while j<=l
    
    temp = Read_Time_Cor(V,U(j), t);
    % Verification that file exist
    if isequal(temp,zeros(1,t))
        
        U(j)=[]; % Does not exist, value for U(j) deleted
        l = l-1; % New length of V array
    else
        % The file exists, 
        plano = cat(2,plano,temp);
        j = j+1;
    end
     
end
%% Plotting parameters
[bU,bT]=meshgrid(U,T);
[DU,Dt] = gradient(plano,deltaH,deltaT); % d/dU & d/dT
L = del2(plano,deltaH,deltaT); % Second derivative
D = diff(plano)/deltaH;
D = cat(1,D,D(end,:));
D2 = diff(D)/deltaH;
D2 = cat(1,D2,D2(end,:));
chosen_times=[2,10];
%% Figure 1: 3D plot of TC vs V and T
figure(1);
surf(bU,bT,plano);
title(['3D Time Correlations for V=' num2str(V)])
ylabel('Tiempo')
xlabel('U')
zlabel('Time Correlations')
%% Figure 2: 3D plot of gradient of TC vs V
figure(2);
surf(bU,bT,plano);
title(['Gradient of TC for V=' num2str(V)])
ylabel('Tiempo')
xlabel('U')
zlabel('$\displaystyle\frac{d}{dU}$','interpreter','latex')
%% Figure 3: 3D plot of Laplacian of TC vs V and T
figure(3)
surf(bU,bT,plano);
title(['Laplacian of TC for V=' num2str(V)])
ylabel('Tiempo')
xlabel('U')
zlabel('\nabla ^2','Interpreter','tex')
%% Figure 4: TC and first 2 derivatives for chosen times
figure(4);
plot(U,plano(chosen_times(1),:),'b','linewidth',2)
hold on
plot(U,plano(chosen_times(2),:),'r','linewidth',2)
plot(U,D(chosen_times(1),:),'g','linewidth',2)
plot(U,D2(chosen_times(1),:),'m','linewidth',2)
title(['Time Correlations for V=' num2str(V) ...
    ' at times ' num2str(T(chosen_times(1))) ' and '...
    num2str(T(chosen_times(2)))])
xlabel('U')
ylabel('Time Correlations (real)')
legend(['T=' num2str(T(chosen_times(1)))],...
    ['T=' num2str(T(chosen_times(2)))],...
    'Primera Derivada', 'Segunda Derivada',...
    'Location', 'Best')
%% Vline and text
vline( 8,'-k')
text(7.75,0.85,'CDW')
text(8.25,0.85,'SDW')
