clear
close all
clc
%% System Parameters

V = -1:0.02:-0.5; U = -2; % 1a
% V = -0.1:0.01:0.1; U = -2; % 1b
% V = -1.6:0.01:-1.4; U = 1; % 2a
% V = -1.3:0.02:-0.7; U = 1; % 2b
% V = 0.2:0.02:0.7; U = 1; % 2c
l = length(V);
deltaH=V(2)-V(1); % V step
deltaT= 0.02; % Timestep
T = 0:deltaT:2; % Time vector
t = length(T); 

%% Reading files to set TC matrix
plano = []; % TC matrix initialization
j=1;
while j<=l
    
    temp = Read_Time_Cor(V(j),U, t);
    % Verification that file exist
    if isequal(temp,zeros(1,t))
        
        V(j)=[]; % Does not exist, value for V(j) deleted
        l = l-1; % New length of V array
    else
        % The file exists, 
        plano = cat(2,plano,temp);
        j = j+1;
    end
     
end
%% Filtering
h=fspecial('average',[7,3]);
h2=[0.5;0.5];
plano = imfilter(plano,h);
plano = imfilter(plano,h2);

%% Plotting parameters
[bV,bT] = meshgrid(V,T); % V and T mesh for surf
[Dv,Dt] = gradient(plano,deltaH,deltaT); % d/dV & d/dT
L = del2(plano,deltaH,deltaT); % Second derivative
% Derivative and 2nd derivative using diff
D = diff(plano)/deltaH;
D = cat(1,D,D(end,:));
D2 = diff(D)/deltaH;
D2 = cat(1,D2,D2(end,:));
chosen_times=[2,101]; % Chosen times for plot 4
%% Figure 1: 3D plot of TC vs V and T
figure(1);
surf(bV,bT,plano);
title(['3D Time Correlations for U=' num2str(U)])
ylabel('Tiempo')
xlabel('V')
zlabel('Time Correlations')
%% Figure 2: 3D plot of gradient of TC vs V
figure(2);
surf(bV,bT,Dv);
title(['Gradient of TC for U=' num2str(U)])
ylabel('Tiempo')
xlabel('V')
zlabel('$\displaystyle\frac{d}{dV}$','interpreter','latex')
%% Figure 3: 3D plot of Laplacian of TC vs V and T
figure(3);
surf(bV,bT,L);
title(['Laplacian of TC for U=' num2str(U)])
ylabel('Tiempo')
xlabel('V')
zlabel('\nabla ^2','Interpreter','tex')
%% Figure 4: TC and first 2 derivatives for chosen times
figure(4);
ax1 = subplot(3,1,1);
plot(V,plano(chosen_times(1),:),'b','linewidth',2)
title(['Time Correlations for U=' num2str(U) ...
    ' at times ' num2str(T(chosen_times(1))) ' and '...
    num2str(T(chosen_times(2)))])
hold on
plot(V,plano(chosen_times(2),:),'r','linewidth',2)
legend(['T=' num2str(T(chosen_times(1)))],...
    ['T=' num2str(T(chosen_times(2)))],'Location', 'Best')
hold off
ax2 = subplot(3,1,2);
plot(V,D(chosen_times(1),:),'g','linewidth',2)
title('First Derivative')
ax3 = subplot(3,1,3);
plot(V,D2(chosen_times(1),:),'m','linewidth',2)
title('Second Derivative')
xlabel('V')
linkaxes([ax1,ax2,ax3],'x')
%% Phase Transitions
% Vline for 1a
vline(-0.75,'-k')
text(-0.9,30,'PS')
text(-0.6,30,'SS') 

% % Vline for 1b
% vline(0,'-k')
% text(-0.05,0.27,'SS')
% text(0.05,0.27,'CDW')

% Vline for 2a
% vline(-1.5,'-k')
% text(-1.55,0.1,'PS')
% text(-1.4,0.1,'TS')

% Vline for 2b
% vline(-0.8,'-k')

% Vline for 2c
% vline([0.4 0.6],{'k-','k-'})
% text(0.3,0.57,'SDW')
% text(0.5,0.57,'BOW')
% text(0.65,0.57,'CDW')