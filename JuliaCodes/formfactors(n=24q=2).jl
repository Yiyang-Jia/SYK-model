include("SYK/unfoldingMethods.jl")
include("SYK/qHermite.jl")
include("SYK/formFactorCompute.jl")
using LinearAlgebra
using DelimitedFiles
using Plots
import StatsPlots 
import PyPlot

evsEns4000 = readdlm("n=24+q=2+ens=4000.txt")
#evsEns4000MATLAB = readdlm("n=24q=2eigenens4000fromMATLAB.txt")

n = 24
q = 2
nens = 4000
stadev = sqrt(factorial(q-1)/(2^q * n^(q-1)))
multiPstadev = sqrt( binomial(n,q) * stadev^2 )

evsEns4000realizations =  reshape(evsEns4000,(2^(Int(n/2)-1),nens))


#spectral form factor before unfolding and subtraction of disconnected piece
times = 0:0.1:100
ssf = zeros(length(times))
for k = 1: length(times)
    allphases = exp.(im*times[k]*evsEns4000realizations)
    phaseColumnSum= sum(eachrow(allphases))
  ssf[k] =  norm(phaseColumnSum)^2/nens 
end

open("n=24q=2SSFdata","w") do  ssffile
    writedlm(ssffile,[collect(times) ssf/2^(n/2-1)])
end

ssfplotdata = readdlm("n=24q=2SSFdata")
plot(log.(ssfplotdata[:,1]), log.(ssfplotdata[:,2]),title="N=24 q=2 spectral form factor", xlabel = "log(t)", ylabel = "log(ssf)",label=nothing)

savefig("n24q2SSF.pdf")


#spectral form factor before unfolding but with disconnected piece subtracted
times = 0:0.1:100
ssfconnected = sffconnectedCompute(times, evsEns4000realizations)

open("n=24q=2SSFconnectedData","w") do  ssfconfile
    writedlm(ssfconfile, ssfconnected)
end

ssfconplotdata = readdlm("n=24q=2SSFconnectedData")
plot(log.(ssfconplotdata[:,1]), log.(ssfconplotdata[:,2]),title="N=24 q=2 connected spectral form factor", xlabel = "log(t)", ylabel = "log(ssf)",label=nothing)
savefig("n24q2SSFconnectedLoglog")

plot((ssfconplotdata[:,1]), log.(ssfconplotdata[:,2]),title="N=24 q=2 connected spectral form factor", xlabel = "t", ylabel = "log(ssf)",label=nothing)
savefig("n24q2SSFconnectedLoglinear")

#spectral form factor unfolded and connected
unfoldedE = readdlm("n=24q=2unfoldedLevels(trapezoidMethod)")
unfoldedRealizations =  reshape(unfoldedE,(2^(Int(n/2)-1),nens))

times = 0:0.01:10
ssfconnectedUnfolded = sffConnectedCompute(times, unfoldedRealizations)

open("n=24q=2SSFconnectedUnfoldedData","w") do  ssfconUnffile
    writedlm(ssfconUnffile, ssfconnectedUnfolded)
end

ssfconUnfplotdata = readdlm("n=24q=2SSFconnectedUnfoldedData")

plot(log.(ssfconUnfplotdata[:,1]), log.(ssfconUnfplotdata[:,2]),title="N=24 q=2 connected unfolded SFF", xlabel = "log(t)", ylabel = "log(ssf)",label=nothing)
savefig("n24q2SSFconnectedUnfoldedLoglog")

plot((ssfconUnfplotdata[:,1]), log.(ssfconUnfplotdata[:,2]),title="N=24 q=2 connected spectral form factor", xlabel = "t", ylabel = "log(ssf)",label=nothing)
savefig("n24q2SSFconnectedUnfoldedLoglinear")

plot((ssfconUnfplotdata[:,1]), (ssfconUnfplotdata[:,2]),title="N=24 q=2 connected spectral form factor", xlabel = "t", ylabel = "ssf",label=nothing)
savefig("n24q2SSFconnectedUnfoldedLinear")


#spectral form factor unfolded and connected at very early time
times = 0: 0.0001: 0.1
unfoldedE = readdlm("n=24q=2unfoldedLevels(trapezoidMethod)")
unfoldedRealizations =  reshape(unfoldedE,(2^(Int(n/2)-1),nens))

ssfEarly = sffConnectedCompute(times, unfoldedRealizations)

open("n=24q=2SSFconnectedUnfoldedEarlyTimeData","w") do  ssfEarlyfile
    writedlm(ssfEarlyfile , ssfEarly)
end

ssfEarlyplotdata = readdlm("n=24q=2SSFconnectedUnfoldedEarlyTimeData")


plot((ssfconUnfplotdata[:,1]), (ssfconUnfplotdata[:,2]),title="N=24 q=2 unfoled SSF early time", xlabel = "t", ylabel = "ssf",label=nothing)
savefig("n24q2SSFearlyTimeLinear")