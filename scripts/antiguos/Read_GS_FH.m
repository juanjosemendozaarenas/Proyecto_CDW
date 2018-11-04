function [O2,O3,E0] = Read_GS_FH(U)
% O2 and O3 are structure factors for a given U
% E0 is the ground state energy
% These 3 values are intialized as NaNs to prevent errors if the
% initfile is not found
O2=NaN;
O3=NaN;
E0=NaN;
%% File parameters
L = 64;
J = 1;
V = 6;
% U = 0;
chi_max = 600;
q = linspace(-pi,pi,L);

fname = ['../output/L64_J1_V6_chi600/GS_FH_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V) '_chi' num2str(chi_max) '.mat'];
check = exist(fname);
if check == 2
    load(fname);
    %% Spin-Density Parameter
    O_sdw=zeros(1,L);
    suma = 0;
    for j=1:L % Loop over q
        for m=1:L
            for n=1:L
                out = szsz_all(m,n) - (nup(m) - ndn(m))*(nup(n)-ndn(n));
                out = exp(-1i*q(j)*(m-n))*out;
                suma = suma + out;
            end
        end
        O_sdw(j) = suma/L; 
    end
    
    % Figure 1
    figure1 = figure(1);
    axes1 = axes('Parent',figure1,...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},...
    'XTick',[-3.1416 -1.5708 0 1.5708 3.1416],...
    'FontSize',20,'FontName','LM Roman Slanted 10');
    xlabel('q')
    ylabel('O_{sdw}')
    box(axes1,'on');
    hold(axes1,'on');
    plot(q,O_sdw)
    ftitle1=['Spin Density parameter U=' num2str(U)];
    title(ftitle1)
    figname1=['../Figures/O_sdp_U_' num2str(U)];
    print(figure1, figname1,'-dpng')
    %% Charge Density Operator
    O_cdw=zeros(1,L);
    
    suma = 0;
    for j=1:L % Loop over q
        for m=1:L
            for n=1:L
                out = nn_all(m,n) - (nup(m) + ndn(m))*(nup(n) + ndn(n));
                out = exp(-1i*q(j)*(m-n))*out;
                suma = suma + out;
            end
        end
        O_cdw(j) = suma/L; 
    end
    % Figure 2
    figure2 = figure(2);
    axes2 = axes('Parent',figure2,...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},...
    'XTick',[-3.1416 -1.5708 0 1.5708 3.1416],...
    'FontSize',20,'FontName','LM Roman Slanted 10');
    box(axes2,'on');
    hold(axes2,'on');
    xlabel('q')
    ylabel('O_{cdw}')
    plot(q,O_cdw)
    ftitle2=['Charge Density parameter U=' num2str(U)];
    title(ftitle2)
    figname2=['../Figures/O_cdp_U_' num2str(U)];
    print(figure2, figname2,'-dpng')
    %% Outputs
    E0 = E(end);
    O2 = O_sdw(end);
    O3 = O_cdw(end);
    
else
    sprintf('File not found')
end

