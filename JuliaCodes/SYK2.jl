# using Pkg
# Pkg.add("LinearAlgebra")
# Pkg.add("Distributions")
# Pkg.add("DelimitedFiles")
# Pkg.add("PyCall")
include("/home/yiyang/Desktop/Codes/julia/DiracMat/dirac.jl")
using LinearAlgebra
using Distributions
using DelimitedFiles
using PyCall
@pyimport  numpy.linalg as nl
#using BenchmarkTools

    

n = 24
q = 2
nens = 4000
stadev = sqrt(factorial(q-1)/(2^q * n^(q-1)))
upperRightList = Array{Any}(missing, n);
for w = 1:n
    upperRightList[w]= diracEvenURblock(n,w);
end

open("n=$n+q=$q+ens=$nens.txt","w") do  file
    

    for j=1:nens
        randBlockHami = sparse(zeros( 2^(Int(n/2) -1), 2^(Int(n/2) -1) ));
        for k1 = 1: n-1
            q2Dirac1 = upperRightList[k1];
            for k2 = k1+1: n
                q2Dirac2 = im* q2Dirac1 * adjoint(upperRightList[k2]);
                coupling = rand(Normal(0,stadev))
                randBlockHami = randBlockHami +coupling*q2Dirac2;
        
            end
        end
    
         ev = nl.eigvalsh(collect(randBlockHami))
         writedlm(file, ev)

        if mod(j,50)==0
            println(j)
        end
    end

end

