function [ plano ] = get_planoV( V, U, T, str )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

l = length(V);
t = length(T);
plano = []; % TC matrix initialization
j=1;
while j<=l
    
    fname = ['output/TE' str '/TimeCorrel_FH' str ...
        '_z_L64_J1_U' num2str(U) '_V' num2str(V(j),'%0.2f')...
        '_chi1200_dt0.01_era0.mat'];
    
    check = exist(fname);
    if check == 2
        load(fname);
        tc = real(TimeCorrel);

    else
        tc = zeros(1,t);
        warning = ['File for U=' num2str(U) ' and V='...
            num2str(V(j)) ' not found. Taking zero por time correlation value'];
        sprintf(warning)
    end
    
    % Verification that file exist
    if isequal(tc,zeros(1,t))
        
        V(j)=[]; % Does not exist, value for V(j) deleted
        l = l-1; % New length of V array
    else
        % The file exists, 
        plano = cat(2,plano,tc);
        j = j+1;
    end
     
end

end
