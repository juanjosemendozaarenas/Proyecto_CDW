%% Creating an initialisation file for calculating ground state of a spin-1/2 XXZ lattice

clear; clc;
path(path,'./matlab/tnt_matfiles/'); % Add path for common functions

%% Define system parameters
L = 8; % System size
Jx = 0; % XX exchange
Jy = -0.5; % YY exchange
Jz_array = [3]; % ZZ exchange
Mag_field = 0; % Magnetic field

use_symm = 0; % Decide if symmetries will be used. 1 if yes, 0 if no
qn_tot = 0; % Zero magnetization

chi_ini_rand = 1; % Truncation parameter for the initial random state without symmetries. Has to be very small so random MPS can be created in C code (this is a bug of the TNT library)
chi = 100; % Initial value of chi
chi_max = 500; % Maximal and chi value of the truncation parameter
delta_chi = 100; % Increase of chi each time maximal error is achieved. Set to zero to keep chi constant.

intermediate = 0; % Set to 1 if initial state of DMRG is a state from previous simulation with same parameters (e.g. if previous simulation was killed for some reason, and intermediate states were being saved)
add_to_file = 0; % Number to add to name of files

% FOR SOME REASON I AM HAVING PROBLEMS WITH THE SYMMETRIES HERE!!! DOES NOT
% HAPPEN FOR BOSONS OR FERMIONS. MIGHT NEED TO CHECK LATER!!! 

rand_wf = 0; % 1 to initialise DMRG with random state created in C, 0 to load from initialisation file
prec = 1e-6; % Precision in calculation of energy
i_max = 60; % Maximal number of DMRG iterations
err_max = 1e-7; % Maximal error before increasing chi in each DMRG step

%% Sweep over parameters
for count_file = 1:size(Jz_array,2)
    
    Jz = Jz_array(count_file);
    
    %% Create n-species 
    ns = 1; % One species of spins
    s = 1/2; % Spin of particles	
    [sx,sy,sz,sp,sm] = tntMatPauliOps(s,ns);
    %[sx,sy,sz,sp,sm] = tntMatSpinOps(s,ns);
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
            qnums(loop,:) = 2*diag(sz{loop})'; % So quantum numbers are integer
        end
        
        tntSystem = tntMatCreateBasisOp(sz{1},qnums);
        tntSystem.sysnum = 1;
        tntSystem.symm_type = 1;
        tntSystem.symm_num_qn = ns;
        
    else % ---------- Use no symmetries ----------
        
        qnums = [];
        
        tntSystem = tntMatCreateBasisOp(sz{1},qnums);
        tntSystem.sysnum = 1; % Type of system (bosonic of spin). Not used by the code.
        tntSystem.symm_type = 0; % No symmetries used
        tntSystem.symm_num_qn = 0; % Zero conserved quantities.
        
    end
    
    %% Global parameters used in linear algebra routines
    % These parameters are used while taking SVDs. They are (in order)
    
    % * The tolerance for zeroing matrix values in automatic blocking
    % * The absolute value below which all singular values will be discarded
    % * The ratio between smallest singular value that will be kept and the largest singular value
    tntSystem.zero_tol = 1e-10;
    tntSystem.abs_trunc_tol = -1;
    tntSystem.rel_trunc_tol = 1e-10;
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
    
    %% Initial state, if loaded from initialisation file
    if(rand_wf == 0)
        
        vac = zeros(d,1); vac(1) = 1; % This corresponds to spin up
        wf = cell(1,L);
        
        % Create an antiferromagnetic state |D U D U D U ... D U>
        for site = 1:L
            if(mod(site,2)==1) % Odd site is down.
                wf{site} = sm{1}*vac;
            elseif(mod(site,2)==0) % Even site is up.
                wf{site} = vac;
            end
        end

        wf = tntMatCreateProdMps(wf,qnums);
    end
    
    
    %% Define XXZ Hamiltonian
    osg = tntMatCreateOpArray({sz{1}}); % On-site operators: local interaction and chemical potential
    osparamg = [Mag_field*ones(1,L)]; % Set Magnetic field
    
    nnlg = tntMatCreateOpArray({sx{1},sy{1},sz{1}}); % Left-side nearest-neighbour ground-state operators
    nnrg = tntMatCreateOpArray({sx{1},sy{1},sz{1}}); % Right-side nearest-neighbour ground-state operators
    nnparamg = [Jx*ones(1,L-1); % Nearest-neighbour exchange parameters
                Jy*ones(1,L-1);
                Jz*ones(1,L-1)]; 
    
    %% Expectation values to take
    ExOp.os_operators = tntMatCreateOpArray({}); % On-site expectation values
    ExOp.os_labels = {};
    ExOp.nn_operators = tntMatCreateOpArray({}); % Two-site nearest neighbours
    ExOp.nn_labels = {};   
    ExOp.cs_operators = tntMatCreateOpArray({}); % Centered-site
    ExOp.cs_labels = {};
    ExOp.ap_operators = tntMatCreateOpArray({}); % All-pairs
    ExOp.ap_labels = {};
    
    % Expectation values for pairs of nearest neighbors, i.e. < A_{j} B_{j+1} c_{k} d_{k+1} >
    ExOp_NN_pairs.os_operators = tntMatCreateOpArray({}); % Must be empty
    ExOp_NN_pairs.os_labels = {}; % Must be empty
    ExOp_NN_pairs.nn_operators = tntMatCreateOpArray({}); % Must be empty
    ExOp_NN_pairs.nn_labels = {}; % Must be empty
    ExOp_NN_pairs.cs_operators = tntMatCreateOpArray({}); % Must be empty
    ExOp_NN_pairs.cs_labels = {}; % Must be empty    
    ExOp_NN_pairs.ap_operators = tntMatCreateOpArray...
        ({sx{1},sx{1},sz{1},sz{1},...
        sx{1},sx{1},sy{1},sy{1},...
        sy{1},sy{1},sz{1},sz{1}}); 
    ExOp_NN_pairs.ap_labels = {'sxsxszsz','sxsxsysy','sysyszsz'};  
    
    %% Save information
    
    % Name of file where ground state will be saved
    savefile = ['GS_XXZ_L' num2str(L) '_Jx' num2str(Jx) '_Jy' num2str(Jy) '_Jz' num2str(Jz) '_B' num2str(Mag_field) '_chi' num2str(chi_max) '.mat'];
    
    % Saving current information
    fname = ['./initfiles/initial_GS_XXZ' num2str(count_file+add_to_file) '.mat'];
    save(fname);
    
end