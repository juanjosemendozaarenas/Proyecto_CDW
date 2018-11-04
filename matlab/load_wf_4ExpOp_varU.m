%% Creating an initialisation file for calculating ground state of an extended Fermi-Hubbard lattice

clear; clc;
path(path,'./matlab/tnt_matfiles/'); % Add path for common functions

%% Define system parameters
L = 64; % System size
J = 1; % Hopping
U = -1:0.05:1; % Nearest-neighbor interaction
V = -0.2; % On-site interaction
chi_max = 1500; % Maximal and chi value of the truncation parameter
add_to_file=0;

for cont=1:length(U)
    lname = ['/output/GS4a/GS_FH4a_L' num2str(L) '_J'...
        num2str(J) '_U' num2str(U(cont),'%.2f') '_V'...
        num2str(V) '_chi' num2str(chi_max) '.mat'];
    %% Load previous wavefunction
    check = exist(lname);
    if check == 2
        load(lname,'wf_g');
    else
        sprintf('File not found')
    end

    %% Create fermionic operators
    ns = 2; % We have two species: fermions with spin up and with spin down
    [n,c,cd,Pc,cdP,sz,sz2,dbl] = tntMatFermionOps;
    d = size(sz{1},1);
    %% Expectation values to take
    ExOp.os_operators = tntMatCreateOpArray({dbl{1},sz{1},n{2}}); % Single-site operators
    ExOp.os_labels = {'dbl','sz','nup'};
    ExOp.nn_operators = tntMatCreateOpArray({}); % Two-site nearest-neighbor operators
    ExOp.nn_labels = {};
    ExOp.cs_operators = tntMatCreateOpArray({}); % Two-site central-site operators
    ExOp.cs_labels = {};
    ExOp.ap_operators = tntMatCreateOpArray({}); % Two-site all pairs operators
    ExOp.ap_labels = {};
    
    % Expectation values for pairs of nearest neighbors, i.e. < A_{j} B_{j+1} c_{k} d_{k+1} >
    ExOp_NN_pairs.os_operators = tntMatCreateOpArray({}); % Must be empty
    ExOp_NN_pairs.os_labels = {}; % Must be empty
    ExOp_NN_pairs.nn_operators = tntMatCreateOpArray({}); % Must be empty
    ExOp_NN_pairs.nn_labels = {}; % Must be empty
    ExOp_NN_pairs.cs_operators = tntMatCreateOpArray({}); % Must be empty
    ExOp_NN_pairs.cs_labels = {}; % Must be empty    
    ExOp_NN_pairs.ap_operators = tntMatCreateOpArray({}); 
    ExOp_NN_pairs.ap_labels = {};  

    %% Save information

    % Name of file Expected values will be saved
    savefile = ['GS_FH4a_L' num2str(L) '_J' num2str(J) '_U'...
        num2str(U(cont),'%.2f') '_V' num2str(V) '_chi'...
        num2str(chi_max) '.mat'];

    % Initfile/existing file name
    fname = ['./initfiles/existing_GS_FH4a_4ExpOps'...
        num2str(cont + add_to_file) '.mat'];
    save(fname);
end