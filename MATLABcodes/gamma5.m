function y = gamma5(n)
if mod(n, 2) == 0
   y = sparse(kron(pauli(3), speye(2^(n/2 -1))));
else 
    disp('gamma5 only exists in even dimensions!')
    
end