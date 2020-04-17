n=16; %number of Majoranas
q=4;
nens=1; %size of the ensemble
stadev = sqrt(factorial(q-1)/(2^q * n^(q-1))); % SYK J=1
upperRightList = cell(1,n);
for i = 1:n
    upperRightList{1, i} = diracEvenURblock(n,i);
end


fileID = fopen('n=16eigenOneRealization.txt','w');
for w=1:nens
    randBlockHami = sparse(zeros(2^(n/2 -1)));
    for k1 = 1: n-3
        q4Dirac1 = upperRightList{1,k1};
        for k2 = k1+1: n-2
            q4Dirac2 = q4Dirac1* ctranspose(upperRightList{1,k2});
            for k3 = k2+1:n-1
                 q4Dirac3 = q4Dirac2*upperRightList{1,k3};
                for k4 = k3+1:n
                    q4Dirac = q4Dirac3* ctranspose(upperRightList{1,k4});
            randBlockHami = randBlockHami +normrnd(0,stadev)*q4Dirac;
                end
            end
        end
    end
  
    ev = eig(full(randBlockHami));
    fprintf(fileID,'%.8f\n', ev);
    if mod(w,10)==0
       w
    end
end
fclose(fileID);

