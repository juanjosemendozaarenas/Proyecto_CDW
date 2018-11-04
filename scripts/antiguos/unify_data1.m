warning('off','MATLAB:load:variableNotFound');
U = -2;
% Defining sub-intervals and properties
V1 = -1:0.02:-0.5; l1 = length(V1);
V2 = -0.45:0.05:-0.15; l2 = length(V2);
V3 = -0.1:0.01:0.1; l3 = length(V3);
V4 = 0.15:0.05:1; l4 = length(V4);
% Interval
V = [V1 V2 V3 V4];

for i=1:length(V)
    % Evaluates the simulation name
    if V(i) < -0.45
        sim = '1a';
    elseif V(i) < 0.55
        sim = '1b';
    else
        sim = '1c';
    end
    %% Loads Time Correlations of dbl
    fname1 = ['/output/TC' sim '_dbl/TimeCorrel_FH' sim ...
        '_dbl_L64_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era0.mat'];
    
    check = exist(fname1);
    if check ~=  2
        fname1 = ['/output/TC' sim '_dbl/TimeCorrel_FH' sim ...
        '_dbl_L64_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era1.mat'];
    end
    
    TC_dbl = load(fname1, 'TimeCorrel');
    TC_dbl = TC_dbl.('TimeCorrel');
    err_dbl = load(fname1, 'err');
    err_dbl = err_dbl.('err');
    load(fname1, 'time');
    load(fname1, 'tstep');
    
    %% Loads Time Correlations of sz
    fname2 = ['/output/TE' sim '/TimeCorrel_FH' sim ...
        '_z_L64_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era0.mat'];
    
    check = exist(fname2);
    if check ~=  2
        fname2 = ['/output/TE' sim '/TimeCorrel_FH' sim ...
        '_z_L64_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era1.mat'];
    end
    
    TC_sz = load(fname2, 'TimeCorrel');
    TC_sz = TC_sz.('TimeCorrel');
    err_sz = load(fname2, 'err');
    err_sz = err_sz.('err');
    
    %% Loads Time Correlations of nup
    fname3 = ['/output/TC' sim '_nup/TimeCorrel_FH' sim ...
        '_nup_L64_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era0.mat'];
    
    check = exist(fname3);
    if check ~=  2
        fname3 = ['/output/TC' sim '_nup/TimeCorrel_FH' sim ...
        '_nup_L64_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era1.mat'];
    end
    
    TC_nup = load(fname3, 'TimeCorrel');
    TC_nup = TC_nup.('TimeCorrel');
    err_nup = load(fname3, 'err');
    err_nup = err_nup.('err');
    
    %% Loads Time Correlations of ndn
    fname3 = ['/output/TC' sim '_ndn/TimeCorrel_FH' sim ...
        '_ndn_L64_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era0.mat'];
    
    check = exist(fname3);
    if check ~=  2
        fname3 = ['/output/TC' sim '_ndn/TimeCorrel_FH' sim ...
        '_ndn_L64_J1_U' num2str(U) '_V' num2str(V(i), '%.2f')...
        '_chi1200_dt0.01_era1.mat'];
    end
    
    TC_ndn = load(fname3, 'TimeCorrel');
    TC_ndn = TC_ndn.('TimeCorrel');
    err_ndn = load(fname3, 'err');
    err_ndn = err_ndn.('err');
    
    %% Load chi and E
    fname4 = ['/output/GS' sim '/GS_FH' sim '_L64_J1_U'...
        num2str(U) '_V' num2str(V(i), '%.2f') '_chi1500.mat'];
    load(fname4,'E');
    load(fname4,'chi');
    load(fname4,'nup');
    check = exist('nup');
    if check == 0
        fname4 = ['/output/ExOps' sim '/GS_FH' sim '_L64_J1_U'...
        num2str(U) '_V' num2str(V(i), '%.2f') '_chi1500.mat'];
        load(fname4,'dbl');
        load(fname4,'sz');
        load(fname4,'nup');
    else 
        load(fname4,'dbl');
        load(fname4,'sz');
    end
    %% Save Variables
    sname = ['/home/juan/Documents/Superconductores/Codigo TNT/' ...
        'output/cmpl_data/data_L64_J1_U'...
        num2str(U) '_V' num2str(V(i), '%.2f') '_chi1500.mat'];
    save(sname, 'chi','E','dbl','sz','nup','time','tstep',...
        'err_dbl','err_sz','err_nup','err_ndn','TC_dbl',...
        'TC_sz','TC_nup','TC_ndn')
%     clear chi E dbl sz nup time tstep err_dbl err_sz err_nup...
%         err_ndn TC_dbl TC_sz TC_nup TC_ndn
end    
    