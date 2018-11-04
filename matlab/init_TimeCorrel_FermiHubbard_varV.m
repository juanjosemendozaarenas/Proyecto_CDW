%% Creating an initialisation file for a TNT simulation of a periodically-driven 1D Fermi-Hubbard lattice 

clear; clc;
path(path,'./matlab/tnt_matfiles/'); % Add path for common functions

%% Define simulation parameters
% Constant parameters of the Hamiltonian
L = 64; % Total number of sites
J = 1; % Hopping
U = -2; % List of values of on-site interaction
V_array = [-0.45:0.05:-0.15 -0.1:0.01:0.1 0.15:0.05:0.5]; % Nearest-neighbour Coulomb interaction

J_coupling = -J*diag(eye(L-1))'; % Array of J
U_coupling = U*diag(eye(L-1))'; % Array of V

add_to_file = 1; % Number to add to name of files
use_symm = 1; % Decide if symmetries will be used. 1 if yes, 0 if no

% Truncation parameters
% chi = 200; % Initial value. SHOULD BE THE MAXIMAL REACHED CHI OF THE DMRG CALCULATION.
update_chi = 1; % Set to 1 to increase chi automatically during simulation
delta_chi = 100; % Increase of chi
chi_max = 1200; % Maximal chi for time evolution
chi_max_gs = 1500; % Maximal chi used for ground state calculation
err_max = 1e-7; % Increase chi if error at each time step is larger than this value

% Parameters for time evolution
intermediate = 0; % Set to 1 if having intermediate simulation
era = 0; % Number of the segment of the total simulation; larger than 1 if intermediate = 1
prev_era = era - 1; % Number of the previous era, used in code only if intermediate = 1 and thus era > 1

T = 2; % Maximum time
dt = 0.01; % Time step
tbigstep = 2; % Calculate expectation values each 'tbigstep' time steps
if(intermediate == 0)
    numsteps = T/dt; % Number of time steps of time evolution
elseif(intermediate == 1)
    T_reached = 1; % Time reached so far
    T_remaining = T-T_reached; % Time remaining
    numsteps = T_remaining/dt;
end

% Information of state saving
save_state_int = 1; % Set to 1 to save evolved state at certain times
savestep = 20; % Save evolved state every savestep steps
save_state_end = 0; % Set to 1 to save final state in main output file

%% Create fermionic operators
ns = 2; % We have two species: fermions with spin up and with spin down
[n,c,cd,Pc,cdP,sz,sz2,dbl] = tntMatFermionOps;
d = size(sz{1},1);

%% Defining the physical basis and symmetry information
% Now give the quantum number(s) for each index of the operator.
% If there are $n$ quantum numbers there should be $n$
% rows in the array, and the number of columns must equal the number of
% rows or columns of the basis operator. Send an empt array if qn
% conservation is not required

if(use_symm == 1) % ---------- Use U(1) symmetry ----------
    
    qnums = zeros(ns,d);
    
    for loop = 1:ns
        qnums(loop,:) = diag(n{loop})';
    end
    
    tntSystem = tntMatCreateBasisOp(n{1},qnums);
    tntSystem.sysnum = 2; % Type of system (bosonic of spin). Not used by the code.
    tntSystem.symm_type = 1; % Using U(1) symmetry
    tntSystem.symm_num_qn = ns; % Number of conserved quantities. Equal to ns e.g. if the number of each species is conserved
    
else % ---------- Use no symmetries ----------
    
    qnums = [];
    
    tntSystem = tntMatCreateBasisOp(n{1},qnums);
    tntSystem.sysnum = 2; % Type of system (bosonic of spin). Not used by the code.
    tntSystem.symm_type = 0; % No symmetries used
    tntSystem.symm_num_qn = 0; % Zero conserved quantities.
    
end

%% Global parameters used in linear algebra routines
% These parameters are used while taking SVDs. They are (in order)

% * The tolerance for zeroing matrix values in automatic blocking
% * The absolute value below which all singular values will be discarded
% * The ratio between smallest singular value that will be kept and the largest singular value
% * The total truncation error per TEBD sweep
tntSystem.zero_tol = 1e-9;
tntSystem.abs_trunc_tol = -1;
tntSystem.rel_trunc_tol = 1e-9;
tntSystem.trunc_err_tol = -1;
%%
% Define the function that will be used for calculating the truncation
% error. Choose from 'sumsquares', '1norm', '1normscaled', '2norm',
% '2normscaled'.
tntSystem.trunc_err_func = 'sumsquares';
%%
% Define the type of SVD to use.
tntSystem.svdtype = 1;
%%
% Define the maximum number of iterations for the iterative eigenvalue
% solver. You may want to change this if you get non-convergance errors.
tntSystem.maxeigiter = 300;
%%
% Determine whether reshape re-use is on or off. It is best to have it on
% (|1|) if you have many tensors of a similar size (i.e. $\chi$ is uniform
% throughout the lattice) and off (|0|) if this is not true
tntSystem.reshape_reuse = 1;

%% Operator for which time correlations will be obtained
correl_oper_type = n{1}; % Dicotomic operator
correl_oper = tntMatCreateOpArray({correl_oper_type}); % Put in array of operators, so it is easy to read it in C code
correl_oper_site = ceil(L/2); % Site where operator is to be applied
correl_oper_text = 'ndn'; % Name given to operator, to identify file

%% Expectation values to take
ExOp.os_operators = tntMatCreateOpArray({}); % Single-site operators
ExOp.os_labels = {};
ExOp.nn_operators = tntMatCreateOpArray({}); % Two-site nearest-neighbor operators
ExOp.nn_labels = {};
ExOp.cs_operators = tntMatCreateOpArray({}); % Two-site central-site operators
ExOp.cs_labels = {};
ExOp.ap_operators = tntMatCreateOpArray({}); % Two-site all pairs operators
ExOp.ap_labels = {};

%% Evolution Hamiltonian for each Coulomb coupling
% Note that all the terms must match the dimensions of the basis operator

for count_file = 1:size(V_array,2)
    
    
    %% Define Extended Fermi-Hubbard Hamiltonian
    V = V_array(count_file);
    V_coupling = V*diag(eye(L-1))';
    
    % On-site terms
    ost = tntMatCreateOpArray({n{2}*n{1}}); % On-site operators: n_up*n_down
    osparamt = [U*diag(eye(L))']; % On-site parameters
    
    % Nearest-neighbour terms
    % This corresponds to c^+_dn x c_dn (RL down), c_dn x c^+_dn (LR down),
    % c^+_up x c_up (RL up), and c_up x c^+_up (LR up)
    
    nnlt = tntMatCreateOpArray({cdP{1},Pc{1},cdP{2},Pc{2},n{2}+n{1}});
    nnrt = tntMatCreateOpArray({c{1},cd{1},c{2},cd{2},n{2}+n{1}});
    nnparamt = [J_coupling;
                J_coupling;
                J_coupling;
                J_coupling;
                V_coupling];
    
    %% Save information
   
    % Name of file where initial state (ground state of same parameters) has been saved
    gsfile = ['GS_FH1b_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V,'%.2f') '_chi' num2str(chi_max_gs) '.mat'];
    
    % Name of file where results of time evolution will be saved
    savefile = ['TimeCorrel_FH1b_' correl_oper_text '_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V,'%.2f') '_chi' num2str(chi_max) '_dt' num2str(dt) '_era' num2str(era) '.mat'];
    
    % Name of files where state of previous era was saved, and where state of current era will be saved
    prev_era_state = ['State_TimeCorrel_FH1b_' correl_oper_text '_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V,'%.2f') '_chi' num2str(chi_max) '_dt' num2str(dt) '_era' num2str(prev_era) '.mat'];
    curr_era_state = ['State_TimeCorrel_FH1b_' correl_oper_text '_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V,'%.2f') '_chi' num2str(chi_max) '_dt' num2str(dt) '_era' num2str(era) '.mat'];
    
    % Saving current information
    fname = ['./initfiles/initial_TimeCorrel_FH1b' num2str(count_file+add_to_file) '.mat'];
    save(fname);
    
end