function y = diracOdd(n,i)
if i < n
    y = diracEven(n-1, i);
else
    y = gamma5(n-1);
end

%wrapped in diracMat