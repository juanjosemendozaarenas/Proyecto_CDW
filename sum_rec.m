function [ S ] = sum_rec( T )
%UNTITLED3 Summary of this function goes here
%   Gives a vector with the sums of all previous entries

L = length(T);
S = zeros(1,L);
suma = 0;
for i=1:L
    suma = suma + T(i);
    S(i) = suma;
end

end

