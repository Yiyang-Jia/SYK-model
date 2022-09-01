#studies the overshoot Abi mentioned, using SYK_2 as an example
include("/home/yiyang/Desktop/Codes/julia/SYK/unfoldingMethods.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/qHermite.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/formFactorCompute.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/timeAveraging.jl")
using LinearAlgebra
using DelimitedFiles
using Plots
import StatsPlots 
import PyPlot

n = 24
q = 2
nens = 4000
stadev = sqrt(factorial(q-1)/(2^q * n^(q-1)))
multiPstadev = sqrt( binomial(n,q) * stadev^2 )
eta = suppressionFactor(n, q) #"SYK/qHermite.jl"

unfoldEnergies = readdlm("n=24q=2unfoldedLevels(trapezoidMethod)")
unfoldEnRealizations = reshape(unfoldEnergies,(2^(Int(n/2)-1),nens))

#take a fraction of eigenvalues symmetricallt from the middle
halfLen = Int(length(unfoldEnRealizations[:,1])/2) 

#fraction1 = 1/2
frac = 1/2
nOneSide = Int( round( frac * halfLen ) )
unfoldEnRealizationsTrunc1 = unfoldEnRealizations[ halfLen - nOneSide : halfLen + nOneSide, : ]

times = 0 : 0.001 :1 
sffConUnfTrunc1 = sffConnectedCompute(times, unfoldEnRealizationsTrunc1)
phalf = plot(sffConUnfTrunc1[: , 1], sffConUnfTrunc1[: , 2])

open("n=24q=2SFFconnectedDataHalfEigenVals","w") do  sffconfileTrunc1
    writedlm(sffconfileTrunc1, sffConUnfTrunc1)
end

#fraction2 = 1/3
frac = 1/3
nOneSide = Int( round( frac * halfLen ) )
unfoldEnRealizationsTrunc2 = unfoldEnRealizations[ halfLen - nOneSide : halfLen + nOneSide, : ]

times = 0 : 0.001 :1 
sffConUnfTrunc2 = sffConnectedCompute(times, unfoldEnRealizationsTrunc2)
pOneThird = plot(sffConUnfTrunc2[: , 1], sffConUnfTrunc2[: , 2])

open("n=24q=2SFFconnectedDataOneThirdEigenVals","w") do  sffconfileTrunc2
    writedlm(sffconfileTrunc2, sffConUnfTrunc2)
end

plot(sffConUnfTrunc1[: , 1], sffConUnfTrunc1[: , 2] - sffConUnfTrunc2[: , 2]) 

#fraction3 = 1/4
frac = 1/4
nOneSide = Int( round( frac * halfLen ) )
unfoldEnRealizationsTrunc3 = unfoldEnRealizations[ halfLen - nOneSide : halfLen + nOneSide, : ]

times = 0 : 0.001 :1 
sffConUnfTrunc3 = sffConnectedCompute(times, unfoldEnRealizationsTrunc3)
pOneQuater = plot(sffConUnfTrunc3[: , 1], sffConUnfTrunc3[: , 2])

open("n=24q=2SFFconnectedDataOneQuarterEigenVals","w") do  sffconfileTrunc3
    writedlm(sffconfileTrunc3, sffConUnfTrunc3)
end

plot(sffConUnfTrunc3[: , 1], sffConUnfTrunc3[: , 2] - sffConUnfTrunc1[: , 2]) 

compare1Plot = plot(sffConUnfTrunc2[80:500], (sffConUnfTrunc2[: , 2] - sffConUnfTrunc1[: , 2])[80:500])
compare = plot!(compare1Plot, sffConUnfTrunc3[80:500], (sffConUnfTrunc3[: , 2] - sffConUnfTrunc1[: , 2])[80:500])

halfwindow = 4
decay1 = [sffConUnfTrunc1[100:end , 1] timeAve(halfwindow , sffConUnfTrunc1[100:end , 2])]
#plot(decay1[:,1] , decay1[:, 2])

decay2 = [sffConUnfTrunc2[100:end , 1] timeAve(halfwindow , sffConUnfTrunc2[100:end , 2])]
#plot(decay2[:,1] , decay2[:, 2])

decay3 = [sffConUnfTrunc3[100:end , 1] timeAve(halfwindow , sffConUnfTrunc3[100:end , 2])]
#plot(decay3[:,1] , decay3[:, 2])


plot(decay1[:,1] , [decay1[:, 2]  decay3[:, 2]] , label = ["half of unfolded spectrum" "a quarter of unfolded spectrum"], title = "SYK_2 (N=24)  SFF overshoot region", xlabel = "t")
#savefig("SYK2overshoot.pdf")

plot(log.(decay1[:,1]) , log.([decay1[:, 2]  decay3[:, 2]]) , label = ["half of unfolded spectrum" "a quarter of unfolded spectrum"], title = "SYK_2 (N=24)  SFF overshoot region loglog plot", xlabel = "t")