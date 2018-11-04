clear
close all
clc
%%
V_array = [-4 -3 linspace(-2,1,20) 2 3 4];
l = length(V_array);
O_all = zeros(1,l);
O_singlet = zeros(1,l);
O_triplet = zeros(1,l);

% Read ExpOp for each value of V
for j=1:l
    [O_all(j),O_singlet(j),O_triplet(j)] = Read_GS_FH_SCParam(V_array(j));
    close all
end
%% Plot O_sdw
figure(1)
axes('FontSize',20,'FontName','LM Roman Slanted 10');
box('on');
plot(V_array,O_singlet)
xlabel('V')
ylabel('O_{snn}')
title('U=1')

%% Plot O_cdw
figure(2)
axes('FontSize',20,'FontName','LM Roman Slanted 10');
box('on');
plot(V_array,O_triplet)
xlabel('V')
ylabel('O_{tnn}')
title('U=1')

%% Plot E0
figure(3)
axes('FontSize',20,'FontName','LM Roman Slanted 10');
box('on');
plot(V_array,O_all)
xlabel('V')
ylabel('O_{sos}')
title('U=1')

