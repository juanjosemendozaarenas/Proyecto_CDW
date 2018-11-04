%% Creating an initialisation file for a TNT simulation of a periodically-driven 1D Fermi-Hubbard lattice 

clear; clc;
path(path,'./tnt_matfiles/'); % Add path for common functions

%% Define simulation parameters
% Constant parameters of the Hamiltonian of time evolution
L = 32; % Total number of sites
J = 1; % Hopping
U_array = [4]; % List of values of on-site interaction
V = 0; % Nearest-neighbor interaction

add_to_file = 0; % Number to add to name of files
use_symm = 1; % Decide if symmetries will be used. 1 if yes, 0 if no

% Select type of state at t = 0
GS = 1; % 1 to read ground state from DMRG result, 0 to use simple product initial state
M = 0; % Magnetization of initial state
if(GS == 0)
    STATE = 'CDW'; % AFM for antiferro, CDW for charge density wave
elseif(GS == 1)
    V_gs = 1.9; % Value of V of the ground state taken as initial state
    STATE = ['Vgs' num2str(V_gs)];
end

% Parameters for time evolution
intermediate = 0; % Set to 1 if having intermediate simulation
era = 1; % Number of the segment of the total simulation; larger than 1 if intermediate = 1
prev_era = era - 1; % Number of the previous era, used in code only if intermediate = 1 and thus era > 1

T = 8; % Maximum time
dt = 0.01; % Time step
tbigstep = 5; % Calculate expectation values each 'tbigstep' time steps
if(intermediate == 0)
    numsteps = T/dt; % Number of time steps of time evolution
    time_grid = 0:dt:T; % Time vector
elseif(intermediate == 1)
    T_reached = 2; % Time reached so far
    T_remaining = T-T_reached; % Time remaining
    numsteps = T_remaining/dt;
    time_grid = T_reached:dt:T; % Time vector
end
time_grid_interp = time_grid(1:end-1) + diff(time_grid)/2; % Intermediate times

% Time-dependence of hopping J(t) = J0*exp(iAcos(omega*t))
A = 1; % Oscillation amplitude
omega = 1; % Oscillation frequency
J_time = -J*exp(1i*A*cos(omega*time_grid_interp)); % Time-dependent hopping at intermediate times

% Truncation parameters
chi = 200; % Initial value
update_chi = 1; % Set to 1 to increase chi automatically during simulation
delta_chi = 10; % Increase of chi
chi_max = 1200; % Maximal chi
err_max = 5e-8; % Increase chi if error at each time step is larger than this value

% Information of state saving
save_state_int = 1; % Set to 1 to save evolved state at certain times
savestep = 100; % Save evolved state every savestep steps
save_state_end = 1; % Set to 1 to save final state in main output file

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
tntSystem.zero_tol = 1e-10;
tntSystem.abs_trunc_tol = -1;
tntSystem.rel_trunc_tol = 1e-10;
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

%% Starting state
if(GS == 0)
    
    vac = zeros(d,1); vac(1) = 1; % Vacuum state
    
    if(STATE == 'AFM') % Create Neel state
        
        % ---------------------------- Neel AF state ----------------------------
        % |D U D U D U ... D U>, so odd sites are down and even sites are up
        for site = 1:L
            if(mod(site,2)==1) % Odd site is down.
                wf_g{site} = cd{1}*vac;
            elseif(mod(site,2)==0) % Even site is up.
                wf_g{site} = cd{2}*vac;
            end
        end
        
    elseif(STATE == 'CDW') % Create charge density wave state
        
        % |0 2 0 2 0 2 ... 0 2>, so odd sites have holes and even sites doublons
        for site = 1:L
            if(mod(site,2)==1) % Odd site is empty
                wf_g{site} = vac;
            elseif(mod(site,2)==0) % Even site has a doublon
                wf_g{site} = cd{2}*cd{1}*vac;
            end
        end
    end
    
    % Save in TNT structure
    wf_g = tntMatCreateProdMps(wf_g,qnums);
end

%% Expectation values to take
ExOp.os_operators = tntMatCreateOpArray({n{1},n{2},dbl{1}}); % Single-site operators
ExOp.os_labels = {'ndn','nup','dbl'};
ExOp.nn_operators = tntMatCreateOpArray({}); % Two-site nearest-neighbor operators
ExOp.nn_labels = {};
ExOp.cs_operators = tntMatCreateOpArray({}); % Two-site central-site operators
ExOp.cs_labels = {};
ExOp.ap_operators = tntMatCreateOpArray({sz{1},sz{1},n{1}+n{2},n{1}+n{2},cd{2}*cd{1},c{1}*c{2}}); % Two-site all pairs operators
ExOp.ap_labels = {'szsz_all','nn_all','pairing_all'};

%% Evolution Hamiltonian for each Coulomb coupling
% Note that all the terms must match the dimensions of the basis operator

for count_file = 1:size(U_array,2)
    
    
    %% Define Extended Fermi-Hubbard Hamiltonian
    U = U_array(count_file);
    
    % On-site terms
    ost = tntMatCreateOpArray({n{2}*n{1}}); % On-site operators: n_up*n_down
    osparamt = [U*diag(eye(L))']; % On-site parameters
    
    % Nearest-neighbour terms
    % This corresponds to c^+_dn x c_dn (RL down), c_dn x c^+_dn (LR down),
    % c^+_up x c_up (RL up), and c_up x c^+_up (LR up)
    
    nnlt = tntMatCreateOpArray({cdP{1},Pc{1},cdP{2},Pc{2},n{1}+n{2}});
    nnrt = tntMatCreateOpArray({c{1},cd{1},c{2},cd{2},n{1}+n{2}});
    
    % This array, containing the constant part of the hopping, is defined
    % here for convenience. Its values are replaced in the main C function
    nnparamt = [-J*diag(eye(L-1))';
        -J*diag(eye(L-1))';
        -J*diag(eye(L-1))';
        -J*diag(eye(L-1))';
        V*diag(eye(L-1))'];
    
    Num_nn_terms = size(nnparamt,1);
    
    %% Save information
    
    if(GS == 1) % If starting from a ground state
        gsfile = ['GS_FH_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V_gs) '_chi' num2str(chi_max) '_M' num2str(M) '.mat']; % Name of file where initial state has been saved
    end
    
    % Name of file where results of time evolution will be saved. 
    savefile = ['TEBD_FH_' STATE '_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V) '_A' num2str(A) '_omega' num2str(omega) '_chi' num2str(chi_max) '_dt' num2str(dt) '_era' num2str(era) '.mat'];
    
    % Name of files where state of previous era was saved, and where state of current era will be saved
    prev_era_state = ['State_TEBD_FH_' STATE '_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V) '_A' num2str(A) '_omega' num2str(omega) '_chi' num2str(chi_max) '_dt' num2str(dt) '_era' num2str(prev_era) '.mat'];
    curr_era_state = ['State_TEBD_FH_' STATE '_L' num2str(L) '_J' num2str(J) '_U' num2str(U) '_V' num2str(V) '_A' num2str(A) '_omega' num2str(omega) '_chi' num2str(chi_max) '_dt' num2str(dt) '_era' num2str(era) '.mat'];
    
    % Saving current information
    fname = ['../initfiles/initial_TEBD_FH' num2str(count_file+add_to_file) '.mat'];
    save(fname);
    
end