function [P, F, L, chi_max, sentinel] = unify_data( sim )
% UNIFY_DATA reads data from files generated from DMRG evolution and time
% correlations and joins them together in a unique file. The data is from
% Fermi-Hubbard hamiltonian, with different values of V and U.
%   'sim' is the number of the simulation. Can be 1, 2, 3,96 or 4. They
%   determine the ubication of the state in the FH phase diagram.
%   The variable 'sentinel' determines the orientation of the data set. If
%   is positive is vertical, and if is negative is horizontal.
% OUTPUTS: P is the varying parameter o the simulation, F is the fixed
% condition. L is the size of the system. Also returns sentinel.
warning('off','MATLAB:load:variableNotFound');
%% Determining the diferent sets of simulations
switch sim
    case 1
        U = -2;
        V = [-1:0.02:-0.5 -0.45:0.05:-0.15 -0.1:0.01:0.1 0.15:0.05:1];
        L = 64;
        chi_max = 1500;
        sentinel = 1;
    case 2
        U = 1;
        V = [-2:0.02:-1.62 -1.6:0.01:-1.31 -1.3:0.02:-0.4 -0.37:0.03:0.17 0.2:0.02:0.9];
        L = 64;
        chi_max = 1500;
        sentinel = 1;
    case 3
        V = 4;
        U = 7.5:0.05:8.5;
        L = 64;
        chi_max = 1000;
        sentinel = -1;
    case 96
        V = 4;
        U = 7.5:0.05:8.5;
        L = 96;
        chi_max = 1000;
        sentinel = -1;
        sim = 3;
    case 4
        V = -0.2;
        U = -1:0.05:1;
        L = 64;
        chi_max = 1500;
        sentinel = -1;
end
%% ----------------------- VERTICAL DATASET -----------------------------%%
if sentinel > 0
    for i=1:length(V)
        %% TC for sz
        fname_TCsz = ['/output/TC_sz' num2str(sim) '/TimeCorrel_FH_z_L'...
            num2str(L) '_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era0.mat'];
    
        TC_sz = load(fname_TCsz, 'TimeCorrel');
        TC_sz = TC_sz.('TimeCorrel');
        err_sz = load(fname_TCsz, 'err');
        err_sz = err_sz.('err');
        %% TC for dbl
        fname_TCdbl = ['/output/TC_dbl' num2str(sim) '/TimeCorrel_FH_dbl_L'...
            num2str(L) '_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
            '_chi1200_dt0.01_era0.mat'];
    
        TC_dbl = load(fname_TCdbl, 'TimeCorrel');
        TC_dbl = TC_dbl.('TimeCorrel');
        err_dbl = load(fname_TCdbl, 'err');
        err_dbl = err_dbl.('err');
        load(fname_TCdbl, 'time');
        load(fname_TCdbl, 'tstep');
        %% TC for nup
        fname_TCnup = ['/output/TC_nup' num2str(sim) '/TimeCorrel_FH_nup_L'...
            num2str(L) '_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
            '_chi1200_dt0.01_era0.mat'];
    
        TC_nup = load(fname_TCnup, 'TimeCorrel');
        TC_nup = TC_nup.('TimeCorrel');
        err_nup = load(fname_TCnup, 'err');
        err_nup = err_nup.('err');
        %% TC fpor ndn
        fname_TCndn = ['/output/TC_ndn' num2str(sim) '/TimeCorrel_FH_ndn_L'...
            num2str(L) '_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
            '_chi1200_dt0.01_era0.mat'];
    
        TC_ndn = load(fname_TCndn, 'TimeCorrel');
        TC_ndn = TC_ndn.('TimeCorrel');
        err_ndn = load(fname_TCndn, 'err');
        err_ndn = err_ndn.('err');
        %% E, chi and Exp Ops
        fname_GS = ['/output/GS' num2str(sim) '/GS_FH_L' num2str(L)...
            '_J1_U' num2str(U) '_V' num2str(V(i), '%.2f') '_chi'...
            num2str(chi_max) '.mat'];
    
        load(fname_GS,'E');
        load(fname_GS,'chi');
        load(fname_GS,'nup');
        check = exist('nup');
        % If expected value not found in GS file, load ExpOp file
        if check == 0
            fname_GS = ['/output/ExOps' num2str(sim) '/GS_FH_L'...
                num2str(L) '_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
                '_chi' num2str(chi_max) '.mat'];
            load(fname_GS,'dbl');
            load(fname_GS,'sz');
            load(fname_GS,'nup');
        else 
            load(fname_GS,'dbl');
            load(fname_GS,'sz');
        end
        %% Save Variables
        sname = ['/home/juan/Documents/Superconductores/Codigo TNT/' ...
            'output/data/data_L' num2str(L) '_J1_U' num2str(U) '_V'...
            num2str(V(i), '%.2f') '_chi' num2str(chi_max) '.mat'];
        save(sname, 'chi','E','dbl','sz','nup','time','tstep',...
            'err_dbl','err_sz','err_nup','err_ndn','TC_dbl',...
            'TC_sz','TC_nup','TC_ndn')
        clear chi E dbl sz nup time tstep err_dbl err_sz err_nup...
            err_ndn TC_dbl TC_sz TC_nup TC_ndn
    end
    % Determines parameter condition and fixed condition
    P = V;
    F = U;
%% --------------------- HORIZONTAL DATASET -----------------------------%%
elseif sentinel < 0
    for i=1:length(U)
        %% TC for sz
        fname_TCsz = ['/output/TC_sz' num2str(sim) '/TimeCorrel_FH_z_L'...
            num2str(L) '_J1_U' num2str(U(i), '%.2f') '_V' num2str(V)...
            '_chi1200_dt0.01_era0.mat'];
    
        TC_sz = load(fname_TCsz, 'TimeCorrel');
        TC_sz = TC_sz.('TimeCorrel');
        err_sz = load(fname_TCsz, 'err');
        err_sz = err_sz.('err');
        %% TC for dbl
        fname_TCdbl = ['/output/TC_dbl' num2str(sim) '/TimeCorrel_FH_dbl_L'...
            num2str(L) '_J1_U' num2str(U(i), '%.2f') '_V' num2str(V)...
            '_chi1200_dt0.01_era0.mat'];
    
        TC_dbl = load(fname_TCdbl, 'TimeCorrel');
        TC_dbl = TC_dbl.('TimeCorrel');
        err_dbl = load(fname_TCdbl, 'err');
        err_dbl = err_dbl.('err');
        load(fname_TCdbl, 'time');
        load(fname_TCdbl, 'tstep');
        %% TC for nup
        fname_TCnup = ['/output/TC_nup' num2str(sim) '/TimeCorrel_FH_nup_L'...
            num2str(L) '_J1_U' num2str(U(i), '%.2f') '_V' num2str(V)...
            '_chi1200_dt0.01_era0.mat'];
    
        TC_nup = load(fname_TCnup, 'TimeCorrel');
        TC_nup = TC_nup.('TimeCorrel');
        err_nup = load(fname_TCnup, 'err');
        err_nup = err_nup.('err');
        %% TC fpor ndn
        fname_TCndn = ['/output/TC_ndn' num2str(sim) '/TimeCorrel_FH_ndn_L'...
            num2str(L) '_J1_U' num2str(U(i), '%.2f') '_V' num2str(V)...
            '_chi1200_dt0.01_era0.mat'];
    
        TC_ndn = load(fname_TCndn, 'TimeCorrel');
        TC_ndn = TC_ndn.('TimeCorrel');
        err_ndn = load(fname_TCndn, 'err');
        err_ndn = err_ndn.('err');
        %% E, chi and Exp Ops
        fname_GS = ['/output/GS' num2str(sim) '/GS_FH_L' num2str(L)...
            '_J1_U' num2str(U(i), '%.2f') '_V' num2str(V) '_chi'...
            num2str(chi_max) '.mat'];
        
        load(fname_GS,'E');
        load(fname_GS,'chi');
        load(fname_GS,'nup');
        check = exist('nup');
        % If expected value not found in GS file, load ExpOp file
        if check == 0
            fname_GS = ['/output/ExOps' num2str(sim) '/GS_FH_L'...
                num2str(L) '_J1_U' num2str(U(i), '%.2f') '_V' num2str(V)...
                '_chi' num2str(chi_max) '.mat'];
            load(fname_GS,'dbl');
            load(fname_GS,'sz');
            load(fname_GS,'nup');
        else 
            load(fname_GS,'dbl');
            load(fname_GS,'sz');
        end
        %% Save Variables
        sname = ['/home/juan/Documents/Superconductores/Codigo TNT/' ...
            'output/data/data_L' num2str(L) '_J1_U' num2str(U(i), '%.2f')...
            '_V' num2str(V) '_chi' num2str(chi_max) '.mat'];
        save(sname, 'chi','E','dbl','sz','nup','time','tstep',...
            'err_dbl','err_sz','err_nup','err_ndn','TC_dbl',...
            'TC_sz','TC_nup','TC_ndn')
        clear chi E dbl sz nup time tstep err_dbl err_sz err_nup...
            err_ndn TC_dbl TC_sz TC_nup TC_ndn
    end
    % Determines parameter condition and fixed condition
    P = U;
    F = V;
end
end

