n=24; %number of Majoranas
q=2;
nens=2000; %size of the ensemble
stadev = sqrt(factorial(q-1)/(2^q * n^(q-1))); % SYK J=1
upperRightList = cell(1,n);
for w = 1:n
    upperRightList{1, w} = diracEvenURblock(n,w);
end

fileID = fopen('n=24q=2eigen.txt','w');
for j=1:nens
    randBlockHami = sparse(zeros(2^(n/2 -1)));
    for k1 = 1: n-1
        q2Dirac1 = upperRightList{1,k1};
        for k2 = k1+1: n
            q2Dirac2 = i*q2Dirac1* ctranspose(upperRightList{1,k2});
           
            randBlockHami = randBlockHami +normrnd(0,stadev)*q2Dirac2;
      
        end
    end
  
    ev = eig(full(randBlockHami));
   fprintf(fileID,'%.8f\n', ev);
    if mod(j,50)==0
        j
    end
end
fclose(fileID);
