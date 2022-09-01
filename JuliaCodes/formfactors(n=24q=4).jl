include("SYK/unfoldingMethods.jl")
include("SYK/qHermite.jl")
include("SYK/formFactorCompute.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/RMTfunctions.jl")
using LinearAlgebra
using DelimitedFiles
using Plots
import StatsPlots 
import PyPlot

evsEns= readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=4Ens=4000FromMATLAB.txt")

n = 24
q = 4
nens = 4000
stadev = sqrt(factorial(q-1)/(2^q * n^(q-1)))
multiPstadev = sqrt( binomial(n,q) * stadev^2 )

evsEns4000realizations =  reshape(evsEns4000,(2^(Int(n/2)-1),nens))


#spectral form factor before unfolding and subtraction of disconnected piece
times = 0:0.1:100
sff = zeros(length(times))
for k = 1: length(times)
    allphases = exp.(im*times[k]*evsEns4000realizations)
    phaseColumnSum= sum(eachrow(allphases))
  sff[k] =  norm(phaseColumnSum)^2/nens 
end

open("n=24q=4SFFdata","w") do  sfffile
    writedlm(sfffile,[collect(times) sff/2^(n/2-1)])
end

sffplotdata = readdlm("n=24q=4SFFdata")
plot(log.(sffplotdata[:,1]), log.(sffplotdata[:,2]),title="N=24 q=4 spectral form factor", xlabel = "log(t)", ylabel = "log(sff)",label=nothing)

savefig("n24q4SFF.pdf")


#spectral form factor before unfolding but with disconnected piece subtracted
times = 0:0.1:100
sffconnected = sffconnectedCompute(times, evsEns4000realizations)

open("n=24q=4SFFconnectedData","w") do  sffconfile
    writedlm(sffconfile, sffconnected)
end

sffconplotdata = readdlm("n=24q=4SFFconnectedData")
plot(log.(sffconplotdata[:,1]), log.(sffconplotdata[:,2]),title="N=24 q=4 connected spectral form factor", xlabel = "log(t)", ylabel = "log(sff)",label=nothing)
savefig("n24q4SFFconnectedLoglog")

plot((sffconplotdata[:,1]), log.(sffconplotdata[:,2]),title="N=24 q=4 connected spectral form factor", xlabel = "t", ylabel = "log(sff)",label=nothing)
savefig("n24q4SFFconnectedLoglinear")

#spectral form factor unfolded and connected
begin
    unfoldedE = readdlm("n=24q=4unfoldedLevels(trapezoidMethod)")
    unfoldedRealizations =  reshape(unfoldedE,(2^(Int(n/2)-1),nens))

    times = 0:0.01:10
    sffconnectedUnfolded = sffConnectedCompute(times, unfoldedRealizations)

    open("n=24q=4SFFconnectedUnfoldedData","w") do  sffconUnffile
        writedlm(sffconUnffile, sffconnectedUnfolded)
    end
end


sffconUnfplotdata = readdlm("n=24q=4SFFconnectedUnfoldedData")

plot(log.(sffconUnfplotdata[:,1]), log.(sffconUnfplotdata[:,2]),title="N=24 q=4 connected unfolded SFF", xlabel = "log(t)", ylabel = "log(sff)",label=nothing)
savefig("n24q4SFFconnectedUnfoldedLoglog")

plot((sffconUnfplotdata[:,1]), log.(sffconUnfplotdata[:,2]),title="N=24 q=4 connected spectral form factor", xlabel = "t", ylabel = "log(sff)",label=nothing)
savefig("n24q4SFFconnectedUnfoldedLoglinear")


goeplot = plot(0:0.01:10, sffGOE, label = "GOE", linewidth =2.5,legend=:bottomright)
plot!(goeplot, (sffconUnfplotdata[:,1]), (sffconUnfplotdata[:,2]),title="Connected spectral form factor", xlabel = "t", ylabel = "sff",label="N=24, q=4 SYK")
savefig("n24q4SFFconnectedUnfoldedLinear.pdf")


#spectral form factor unfolded and connected at very early time
times = 0: 0.0001: 0.1
unfoldedE = readdlm("n=24q=4unfoldedLevels(trapezoidMethod)")
unfoldedRealizations =  reshape(unfoldedE,(2^(Int(n/2)-1),nens))

sffEarly = sffConnectedCompute(times, unfoldedRealizations)

open("n=24q=4SFFconnectedUnfoldedEarlyTimeData","w") do  sffEarlyfile
    writedlm(sffEarlyfile , sffEarly)
end

sffEarlyplotdata = readdlm("n=24q=4SFFconnectedUnfoldedEarlyTimeData")


plot((sffconUnfplotdata[:,1]), (sffconUnfplotdata[:,2]),title="N=24 q=4 unfoled SFF early time", xlabel = "t", ylabel = "sff",label=nothing)
savefig("n24q42SFFearlyTimeLinear")


#spectral form factor unfolded and connected, 6th order ensemble unfolding
begin
    unfoldedE = readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=4unfoldedLevels(6thOrderEnsembleUnfoldingWithOddTerms)0p9Bulk")
    unfoldedRealizations =  reshape(unfoldedE,(2^(Int(n/2)-1),nens))

    times = 0:0.01:10
    sffconnectedUnfolded = sffConnectedCompute(times, unfoldedRealizations)

    open("n=24q=4SFFconnectedUnfoldedData6thOrderEnsUnfolding","w") do  sffconUnffile
        writedlm(sffconUnffile, sffconnectedUnfolded)
    end
end
goeplot = plot(0:0.01:10, sffGOE, label = "GOE", linewidth =2.5,legend=:bottomright)
plot!(goeplot, sffconnectedUnfolded[:,1], sffconnectedUnfolded[:,2], label="N=24, q=4, 6th order unfolding",legend=:bottomright)

