clear
close all
clc
%%
U_array = [0 2 4 6 8 10:0.25:14 16 18 20];
l = length(U_array);
E0 = zeros(1,l);
O2 = zeros(1,l);
O3 = zeros(1,l);

for j=1:l
    [O2(j),O3(j),E0(j)] = Read_GS_FH(U_array(j));
    close all
end
%% Plot O_sdw
figure(1)
axes('FontSize',20,'FontName','LM Roman Slanted 10');
box('on');
plot(U_array,O2)
xlabel('U')
ylabel('O_{sdw}')

%% Plot O_cdw
figure(2)
axes('FontSize',20,'FontName','LM Roman Slanted 10');
box('on');
plot(U_array,O3)
xlabel('U')
ylabel('O_{cdw}')

%% Plot E0
figure(3)
axes('FontSize',20,'FontName','LM Roman Slanted 10');
box('on');
plot(U_array,E0)
xlabel('U')
ylabel('E_0')

