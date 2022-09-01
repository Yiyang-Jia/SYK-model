#using Pkg
#Pkg.add("RandomMatrices")
using RandomMatrices # documentation https://github.com/JuliaMath/RandomMatrices.jl
include("/home/yiyang/Desktop/Codes/julia/SYK/unfoldingMethods.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/qHermite.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/formFactorCompute.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/timeAveraging.jl")
using LinearAlgebra
using DelimitedFiles
using Plots
import StatsPlots 
import PyPlot


GOEdistr = Wigner{1}()
l = 2048
nens = 10^5
#generate data:
# open("GOE+L=$l+ens=$nens", "w") do file 
#     for i = 1 : nens
#         ev =  eigvalrand(GOEdistr, l)
#         writedlm(file, ev)
#         if mod(i, 500) ==0
#             println(i)
#         end
#     end
# end


function semiCircileDen(x) 
    if  abs(x) < 1
        return l * 2/pi *sqrt(1 - x^2)
    else
        return 0
    end
end

#generate unfolded levels and write to file
# goeEigs = readdlm("GOE+L=$l+ens=$nens")
# goeEigRealizations = reshape(goeEigs, (l, nens))
# goeEigs = nothing #free memory
# unfoldedEigRealizations =  zeros(l, nens)
# for k = 1 : nens
#     unfoldedEigRealizations[:, k] = unfoldTrapz(semiCircileDen , goeEigRealizations[: , k])
# end
# goeEigRealizations = nothing #free memory
# open("GOEUnfolded+L=$l+ens=$nens","w") do unfoFile
#     writedlm(unfoFile, unfoldedEigRealizations[1: l* nens])
# end

# unfoldedEigRealizationsFile = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/GOEUnfolded+L=2048+ens=100000"), (2048, 100000))
# times = 0: 0.01 : 10
# goeSFF = sffConnectedCompute(times , unfoldedEigRealizationsFile)
# open("GOE+L=2048SFFdata","w") do sffFile
#     writedlm(sffFile,goeSFF)
# end


goeSFF = readdlm("/home/yiyang/Desktop/Codes/julia/GOE+L=2048SFFdata")
plot(goeSFF[:,1], goeSFF[:,2])



#Take a fraction (frac) of eigenvalues symmetrically from the middle of the spectrum, the compute SFF from them
# halfLen = Int(length(unfoldedEigRealizationsFile[:,1])/2) 
# frac = 1/2
# nOneSide = Int( round( frac * halfLen ) )
# unfoldEnRealizationsTrunc1 = unfoldedEigRealizationsFile[ halfLen - nOneSide : halfLen + nOneSide, : ]
# unfoldedEigRealizationsFile = nothing #free memory
# times = 0 : 0.01 :10 
# sffConUnfTrunc1 = sffConnectedCompute(times, unfoldEnRealizationsTrunc1)

# open("GOE+L=2048SFFdataHalfEigenVals","w") do  sffconfileTrunc1
#     writedlm(sffconfileTrunc1, sffConUnfTrunc1)
# end
sffHalfSepc = readdlm("/home/yiyang/Desktop/Codes/julia/GOE+L=2048SFFdataHalfEigenVals")
phalf = plot(sffHalfSepc[: , 1], sffHalfSepc[: , 2])

pCompare = plot(sffHalfSepc[: , 1] , [sffHalfSepc[: , 2] goeSFF[:,2]], label = nothing, xlabel="t", title = "SFF of GOE (L=2048, ens=10^5)")

#result: SFF computed from full spectrum and half spectrum of GOE are basically identical