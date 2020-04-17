function y = diracEven(n, i)%Dirac matrices in Weyl basis
if n < 4
    disp('only works for dimension higher than or equal to 4');
elseif i == 1
    y = sparse(kron(pauli(i),speye(2^(n/2 -1))));
elseif i == n
    y = pauli(2);
    for c = 1 : n/2 -1
        y = sparse(kron(y, pauli(2)));
    end
else
    y = pauli(2);
    for c = 1 : floor(i/2) - 1
        y = sparse(kron(y, pauli(2)));
    end
    y = sparse(kron(y, pauli(remod(i))));
    y = sparse(kron(y,speye(2^(n/2 - floor((i)/2)-1))));
end
end
% wrapped in diracMat