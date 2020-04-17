function y = diracMat(n, i)
if mod(n, 2) == 0
    y = diracEven(n, i);
else
    y= diracOdd(n, i);
end