function y = diracEvenURblock(n,i)
%generate the upper-right block of Dirac matrices in even dimensions
%  Our diracEven is already in chiral basis. 
y = diracEven(n, i);
y = y(1:2^(n/2-1),2^(n/2-1)+1:2^(n/2));
end

