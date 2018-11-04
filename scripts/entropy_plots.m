function [ ] = entropy_plots( sim )
% ENTROPY_PLOTS Lee información producida por el script Entanglement.c
% y bota gráficas de Gap de Schmidt contra parametro y ultima
% entrada de entropia contra parametro
%   sim es el nombre de la simulacion, es un String y puede ser '1b',
%   '1c', '2b', '2c', '2d', '3a', '3a96' o '4a'


switch sim
    case '1b'
        L = 64;
        chi = 1500;
        U = -2;
        V = [0:0.01:0.1 0.15:0.05:0.5];
        l = length(V);
        SG = zeros(1,l);
        E = zeros(1,l);
        for i=1:l
            fname = ['./output/Entang' sim '/Entang_chi500_GS_FH'...
                sim '_L' num2str(L) '_J1_U' num2str(U) '_V'...
                num2str(V(i),'%0.2f') '_chi' num2str(chi) '.mat'];
            load(fname, 'Schmidt_gap','Entropy');
            SG(i) = Schmidt_gap;
            E(i) = Entropy(end);
        end
        figure(1)
        plot(V,SG,'o')
        title('Gap de Schmidt para U = -2')
        xlabel('V')
        ylabel('Schmidt Gap');
        
        figure(2)
        plot(V,E,'o')
        title('Última entrada de entropía para U = -2')
        xlabel('V')
        ylabel('E(end)');  
    case '1c'
        L = 64;
        chi = 1500;
        U = -2;
        V = 0.55:0.05:1;
        l = length(V);
        SG = zeros(1,l);
        E = zeros(1,l);
        for i=1:l
            fname = ['./output/Entang' sim '/Entang_chi500_GS_FH'...
                sim '_L' num2str(L) '_J1_U' num2str(U) '_V'...
                num2str(V(i),'%0.2f') '_chi' num2str(chi) '.mat'];
            load(fname, 'Schmidt_gap','Entropy');
            SG(i) = Schmidt_gap;
            E(i) = Entropy(end);
        end
        figure()
        plot(V,SG,'o')
        title('Gap de Schmidt para U = -2')
        xlabel('V')
        ylabel('Schmidt Gap');
        
        figure()
        plot(V,E,'o')
        title('Última entrada de entropía para U = -2')
        xlabel('V')
        ylabel('E(end)');
    case '2b'
        L = 64;
        chi = 1500;
        U = 1;
        V = -0.48:0.02:-0.4;
        l = length(V);
        SG = zeros(1,l);
        E = zeros(1,l);
        for i=1:l
            fname = ['./output/Entang' sim '/Entang_chi500_GS_FH'...
                sim '_L' num2str(L) '_J1_U' num2str(U) '_V'...
                num2str(V(i),'%0.2f') '_chi' num2str(chi) '.mat'];
            load(fname, 'Schmidt_gap','Entropy');
            SG(i) = Schmidt_gap;
            E(i) = Entropy(end);
        end
        figure(1)
        plot(V,SG,'o')
        title('Gap de Schmidt para U = 1')
        xlabel('V')
        ylabel('Schmidt Gap');
        
        figure(2)
        plot(V,E,'o')
        title('Última entrada de entropía para U = 1')
        xlabel('V')
        ylabel('E(end)');
    case '2c'
        L = 64;
        chi = 1500;
        U = 1;
        V = 0.2:0.02:0.9;
        l = length(V);
        SG = zeros(1,l);
        E = zeros(1,l);
        for i=1:l
            fname = ['./output/Entang' sim '/Entang_chi500_GS_FH'...
                sim '_L' num2str(L) '_J1_U' num2str(U) '_V'...
                num2str(V(i),'%0.2f') '_chi' num2str(chi) '.mat'];
            load(fname, 'Schmidt_gap','Entropy');
            SG(i) = Schmidt_gap;
            E(i) = Entropy(end);
        end
        figure(1)
        plot(V,SG,'o')
        title('Gap de Schmidt para U = 1')
        xlabel('V')
        ylabel('Schmidt Gap');
        
        figure(2)
        plot(V,E,'o')
        title('Última entrada de entropía para U = 1')
        xlabel('V')
        ylabel('E(end)');
    case '2d'
        L = 64;
        chi = 1500;
        U = 1;
        V = -0.37:0.03:0.17;
        l = length(V);
        SG = zeros(1,l);
        E = zeros(1,l);
        for i=1:l
            fname = ['./output/Entang' sim '/Entang_chi500_GS_FH'...
                sim '_L' num2str(L) '_J1_U' num2str(U) '_V'...
                num2str(V(i),'%0.2f') '_chi' num2str(chi) '.mat'];
            load(fname, 'Schmidt_gap','Entropy');
            SG(i) = Schmidt_gap;
            E(i) = Entropy(end);
        end
        figure(1)
        plot(V,SG,'o')
        title('Gap de Schmidt para U = 1')
        xlabel('V')
        ylabel('Schmidt Gap');
        
        figure(2)
        plot(V,E,'o')
        title('Última entrada de entropía para U = 1')
        xlabel('V')
        ylabel('E(end)');
    case '3a'
        L = 64;
        chi = 1000;
        V = 4;
        U = 7.5:0.05:8.5;
        l = length(U);
        SG = zeros(1,l);
        E = zeros(1,l);
        for i=1:l
            fname = ['./output/Entang' sim '/Entang_chi500_GS_FH'...
                sim '_L' num2str(L) '_J1_U' num2str(U(i),'%0.2f')...
                '_V' num2str(V) '_chi' num2str(chi) '.mat'];
            load(fname, 'Schmidt_gap','Entropy');
            SG(i) = Schmidt_gap;
            E(i) = Entropy(end);
        end
        figure(1)
        hold on
        plot(U,SG,'o')
        title('Gap de Schmidt para V = 4')
        xlabel('U')
        ylabel('Schmidt Gap');
        
        figure(2)
        hold on
        plot(U,E,'o')
        title('Última entrada de entropía para V = 4')
        xlabel('U')
        ylabel('E(end)');
    case '3a96'
        L = 96;
        chi = 1000;
        V = 4;
        U = 7.5:0.05:8.5;
        l = length(U);
        SG = zeros(1,l);
        E = zeros(1,l);
        for i=1:l
            fname = ['./output/Entang' sim '/Entang_chi500_GS_FH3a'...
                '_L' num2str(L) '_J1_U' num2str(U(i),'%0.2f')...
                '_V' num2str(V) '_chi' num2str(chi) '.mat'];
            load(fname, 'Schmidt_gap','Entropy');
            SG(i) = Schmidt_gap;
            E(i) = Entropy(end);
        end
        figure(1)
        plot(U,SG,'o')
        title('Gap de Schmidt para V = 4, L=96')
        xlabel('U')
        ylabel('Schmidt Gap');
        
        figure(2)
        plot(U,E,'o')
        title('Última entrada de entropía para V = 4, L=96')
        xlabel('U')
        ylabel('E(end)');
    case '4a'
        L = 64;
        chi = 1500;
        V = -0.2;
        U = 0.45:0.05:1;
        l = length(U);
        SG = zeros(1,l);
        E = zeros(1,l);
        for i=1:l
            fname = ['./output/Entang' sim '/Entang_chi500_GS_FH'...
                sim '_L' num2str(L) '_J1_U' num2str(U(i),'%0.2f')...
                '_V' num2str(V) '_chi' num2str(chi) '.mat'];
            load(fname, 'Schmidt_gap','Entropy');
            SG(i) = Schmidt_gap;
            E(i) = Entropy(end);
        end
        figure(1)
        plot(U,SG,'o')
        title('Gap de Schmidt para V = -0.2')
        xlabel('U')
        ylabel('Schmidt Gap');
        
        figure(2)
        plot(U,E,'o')
        title('Última entrada de entropía para V = -0.2')
        xlabel('U')
        ylabel('E(end)');
    otherwise
        msg = ['Error: ' sim 'is not a valid input.'];
        error(msg)
end


    


end

