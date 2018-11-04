%% Creating an initialisation file for calculating ground state of an extended Fermi-Hubbard lattice

clear; clc;
path(path,'./matlab/tnt_matfiles/'); % Add path for common functions

%% Define system parameters
L = 64; % System size
J = 1; % Hopping
U_array = -1:0.05:1; % List of values of on-site interaction
V = -0.2; % Nearest-neighbor interaction

J_coupling = -J*diag(eye(L-1))'; % Array of J
V_coupling = V*diag(eye(L-1))'; % Array of V

% Mag=0 % Give an initial magnetization
% qn_tot = [0.5*L-Mag 0.5*L+Mag]; % Define quantum numbers for each species
qn_tot = [0.5*L 0.5*L]; % Define quantum numbers for each species
add_to_file = 0; % Number to add to name of files

chi_ini_rand = 1; % Truncation parameter for the initial random state without symmetries. Has to be very small so random MPS can be created in C code (this is a bug of the TNT library)
chi = 100; % Initial value of chi
chi_max = 1500; % Maximal and chi value of the truncation parameter
delta_chi = 100; % Increase of chi each time maximal error is achieved. Set to zero to keep chi constant.

intermediate = 0; % Set to 1 if initial state of DMRG is a state from previous simulation with same parameters (e.g. if previous simulation was killed for some reason, and intermediate states were being saved)
use_symm = 1; % Decide if symmetries will be used. 1 if yes, 0 if no

rand_wf = 1; % 1 to initialise DMRG with random state created in C, 0 to load from initialisation file
prec = 1e-6; % Precision in calculation of energy
i_max = 60; % Maximal number of DMRG iterations
err_max = 1e-7; % Maximal error before increasing chi in each DMRG step

%% Sweep over parameters
for count_file = 1:size(U_array,2)
    
    U = U_array(count_file);
    
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
            qnums(loop,:) = diag(n{loop})'; % Quantum numbers to be used per site
        end
        
        tntSystem = tntMatCreateBasisOp(n{1},qnums);
        tntSystem.sysnum = 2; % Type of system (bosonic of spin). Not used by the code.
        tntSystem.symm_type = 1; % Using U(1) symmetry
        tntSystem.symm_num_qn = ns; % Number of conserved quantities
        
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
    tntSystem.zero_tol = 1e-9;
    tntSystem.abs_trunc_tol = -1;
    tntSystem.rel_trunc_tol = 1e-9;
    tntSystem.trunc_err_tol = -1;
    
    % Define the function that will be used for calculating the truncation
    % error. Choose from 'sumsquares', '1norm', '1normscaled', '2norm',
    % '2normscaled'.
    tntSystem.trunc_err_func = 'sumsquares';
    
    % Define the type of SVD to use.
    tntSystem.svdtype = 1;
    
    % Define the maximum number of iterations for the iterative eigenvalue
    % solver. You may want to change this if you get non-convergance errors.
    tntSystem.maxeigiter = 300;
    
    % Determine whether reshape re-use is on or off. It is best to have it on
    % (|1|) if you have many tensors of a similar size (i.e. $\chi$ is uniform
    % throughout the lattice) and off (|0|) if this is not true
    tntSystem.reshape_reuse = 1;
    
    %% Define Extended Fermi-Hubbard Hamiltonian
    U_coupling = U*diag(eye(L))';
    
    % On-site terms
    osg = tntMatCreateOpArray({n{2}*n{1}}); % On-site operators: n_up*n_down
    osparamg = [U_coupling]; % On-site parameters
    
    % Nearest-neighbour terms
    % This corresponds to c^+_dn x c_dn (RL down), c_dn x c^+_dn (LR down),
    % c^+_up x c_up (RL up), c_up x c^+_up (LR up), and n x n (Nearest-neighbour coupling)
    
    nnlg = tntMatCreateOpArray({cdP{1},Pc{1},cdP{2},Pc{2},n{1}+n{2}});
    nnrg = tntMatCreateOpArray({c{1},cd{1},c{2},cd{2},n{1}+n{2}});

    % This array, containing the constant part of the hopping, is defined
    % here for convenience. Its values are replaced in the main C function
    nnparamg = [J_coupling;
        J_coupling;
        J_coupling;
        J_coupling;
        V_coupling];
    
    %% Expectation values to take
    ExOp.os_operators = tntMatCreateOpArray({}); % Single-site operators
    ExOp.os_labels = {};
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
    
    %% Define initial state of DMRG simulation, if it is not created randomly in the C code
    % The filling of this state defines the filling of the ground state
    
    if(rand_wf==0 && intermediate==0) % Initial AF product state (half-filling)
        
        Filling = 'Half'; % Select type of filling        
        vac = zeros(d,1); vac(1) = 1; % Vacuum state
        
        if (Filling == 'Half') % Initial state with half filling
            
            % |D U D U D U ... D U>, AF state so odd sites are down and even sites are up
            for site = 1:L
                if(mod(site,2)==1) % Odd site is down.
                    wf{site} = cd{1}*vac;
                elseif(mod(site,2)==0) % Even site is up.
                    wf{site} = cd{2}*vac;
                end
            end       
        end
        
        % Save in TNT structure
        wf = tntMatCreateProdMps(wf,qnums);        
    end
    
    %% Save information
    
    % Name of file where ground state will be saved
    savefile = ['GS_FH4a_L' num2str(L) '_J' num2str(J) '_U' num2str(U,'%0.2f') '_V' num2str(V) '_chi' num2str(chi_max) '.mat'];
    
    % Saving current information
    fname = ['./initfiles/initial_GS_FH4a' num2str(count_file+add_to_file) '.mat'];
    save(fname);
    
end