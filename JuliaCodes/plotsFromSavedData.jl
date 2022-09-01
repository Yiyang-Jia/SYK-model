begin
    include("/home/yiyang/Desktop/Codes/julia/SYK/numberVarianceCompute.jl")
    include("/home/yiyang/Desktop/Codes/julia/SYK/RMTfunctions.jl")
    include("/home/yiyang/Desktop/Codes/julia/SYK/unfoldingMethods.jl")
    using DelimitedFiles
    using Plots
    import StatsPlots 
    import PyPlot
    using LaTeXStrings
    using StatsBase
end

numberVarq2LocalUnfolding = readdlm("/home/yiyang/Desktop/Codes/julia/numberVarianceN24q2SixthOrderLocalUnfolding")
numberVarq2EnsembleUnfolding = readdlm("/home/yiyang/Desktop/Codes/julia/numberVarianceN24q2SixthOrderEnsembleUnfolding")
nVxdata = numberVarq2LocalUnfolding[:,1]
nVydataLocalUnf = numberVarq2LocalUnfolding[:,2]
nVydataEnsembleUnf = numberVarq2EnsembleUnfolding[:,2]
begin
firstNpoints = 10
plot(nVxdata[1:firstNpoints], [nVydataLocalUnf[1:firstNpoints], nVydataEnsembleUnf[1:firstNpoints],nVxdata[1:firstNpoints]], label = ["6th order local unfoldig" "6th order ensemble unfodling" "Poisson"], leg =:bottomright)
end