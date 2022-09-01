include("/home/yiyang/Desktop/Codes/julia/SYK/unfoldingMethods.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/qHermite.jl")
using LinearAlgebra
using DelimitedFiles
using Plots
import StatsPlots 
import PyPlot

#q = 2
evsq4Ens4000 = readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=4Ens=4000FromMATLAB.txt")
#evsEns4000MATLAB = readdlm("n=24q=2eigenens4000fromMATLAB.txt")

n = 24
q = 4
nens = 4000
stadev = sqrt(factorial(q-1)/(2^q * n^(q-1)))
multiPstadev = sqrt( binomial(n,q) * stadev^2 )
eta = suppressionFactor(n, q) #"SYK/qHermite.jl"

#specDenq2(x) = 2^(n/2-1)*exp(-x^2/(2*multiPstadev^2))/sqrt(2 * π * multiPstadev^2) #normalize the density so that the integral = total number of states
specDenq4(x) = 2^(n/2-1) * qHdensityNormalized(x,eta, multiPstadev , 30)  #qHdensityNormalized(x,eta, sig, ntrunc) 

evsEns4000realizations =  reshape(evsq4Ens4000,(2^(Int(n/2)-1),nens))
#spectral form factor after unfolding and subtraction of disconnected piece
#testunfold = unfoldRect(specDenq2, evsEns4000realizations[:,1])
#StatsPlots.histogram(testunfold)
    # @time begin
    # unfoldedEnergyLevels =  Array{Any}(nothing, size(evsEns4000realizations))
    # for k = 1 : nens
    #     unfoldedEnergyLevels[:,k] =unfoldRect(specDenq2,evsEns4000realizations[:,k])
    # end 
    # end
    # StatsPlots.histogram(unfoldedEnergyLevels[1:2048*4000], label = nothing)
#trapezoid method
unfoldedEnergyLevels =  Array{Float64}(UndefInitializer(), size(evsEns4000realizations))
for k = 1 : nens
    unfoldedEnergyLevels[:,k] =unfoldTrapz(specDenq4,evsEns4000realizations[:,k])
    if mod(k, 200) == 0
        println(k)
    end
end 


open("n=24q=4unfoldedLevels(trapezoidMethod)","w") do  uftrapfile
writedlm(uftrapfile,unfoldedEnergyLevels[1:2048*4000])
end


unfoldedEnergyLevels = readdlm("n=24q=4unfoldedLevels(trapezoidMethod)")
StatsPlots.histogram(unfoldedEnergyLevels[1:2048*4000], nbins=200, label = "200 bins, 4000 realizations", legend=:bottomright, title = "Unfolded spectral density (N=24, q=4)")
savefig("n=24q=4UnfoledHisto(trapezoidMethod).pdf")


#hquadrature method #Not visibly better than trapezoid method
# unfoldedEnergyLevelshquad =  Array{Any}(nothing, size(evsEns4000realizations))
# for k = 1 : nens
#     unfoldedEnergyLevelshquad[:,k] =unfoldhquad(specDenq2,evsEns4000realizations[:,k])
#     if mod(k, 200) == 0
#         println(k)
#     end
# end 


# StatsPlots.histogram(unfoldedEnergyLevelshquad[1:2048*4000], nbins=150, label = nothing, title = "Unfolded spectral density (N=24, q=2), hquadrature method")
# savefig("n=24q=2UnfoledHisto(hquadratureMethod)")

# open("n=24q=2unfoldedLevels(hquadMethod)","w") do  ufhquadfile
# writedlm(ufhquadfile,unfoldedEnergyLevelshquad[1:2048*4000])
# end

