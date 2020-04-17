function y = pauli(x)
if x == 1
    y = [0, 1; 1, 0];
elseif x == 2
    y = [0, -j; j, 0];
elseif x == 3
    y = [1, 0; 0, -1];
else
    disp('there are only 3 pauli matrices!')
end