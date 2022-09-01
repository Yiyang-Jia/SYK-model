include("/home/yiyang/Desktop/Codes/julia/SYK/unfoldingMethods.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/qHermite.jl")
using LinearAlgebra
using DelimitedFiles
using Plots
using StatsBase
using LsqFit
import StatsPlots 
import PyPlot

n = 24
q = 4
nens = 4000
evsq4Ens4000 = readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=4Ens=4000FromMATLAB.txt")
evsEns4000realizations =  reshape(evsq4Ens4000,(2^(Int(n/2)-1),nens))
#evsEns4000MATLAB = readdlm("n=24q=2eigenens4000fromMATLAB.txt")


stadev = sqrt(factorial(q-1)/(2^q * n^(q-1)))
multiPstadev = sqrt( binomial(n,q) * stadev^2 )
eta = suppressionFactor(n, q) #"SYK/qHermite.jl"
emax = groundStateE(multiPstadev, eta)


function sixthOrderFit(n::Int , eigenData::Array, multiPstadev::Real, eta::Real, emax::Real)
    model3(x, p) = p[7] * qHdensityNormalized.(x,eta, multiPstadev , 30) .* .+(1,  p[1]*qHpolynomial.(x, 1, eta, multiPstadev),  p[2]*qHpolynomial.(x, 2, eta, multiPstadev) , p[3]*qHpolynomial.(x, 3, eta, multiPstadev), p[4]*qHpolynomial.(x, 4, eta, multiPstadev) , p[5]*qHpolynomial.(x, 5, eta, multiPstadev),  p[6]*qHpolynomial.(x, 6, eta, multiPstadev) )
    p0model3 = [0, 0, 0, 0,0, 0, 2^(n/2-1)*nens ]
    nbins = 71
    bulkFraction =0.9
    (xdata, ydata) = histogramData(eigenData, -emax*bulkFraction, emax*bulkFraction, nbins)

    histoCurveQHLsqFitModel3 = curve_fit(model3, xdata, ydata, p0model3) #exclude edges of the histogram when fitting
    fittedCoefModel3= coef(histoCurveQHLsqFitModel3)
    
    #chiSquareModel3 = sum( histoCurveQHLsqFitModel3.resid.^2 ) #chi square if the density were normalized to 1
    fittedDen3(x) =   2^(n/2-1) * qHdensityNormalized(x,eta, multiPstadev , 30) * (1 + fittedCoefModel3[1]*qHpolynomial(x, 1, eta, multiPstadev) + fittedCoefModel3[2]*qHpolynomial(x, 2, eta, multiPstadev) + fittedCoefModel3[3]*qHpolynomial(x, 3, eta, multiPstadev) + fittedCoefModel3[4]*qHpolynomial(x, 4, eta, multiPstadev) + fittedCoefModel3[5]*qHpolynomial(x, 5, eta, multiPstadev) +  fittedCoefModel3[6]*qHpolynomial(x, 6, eta, multiPstadev) )
    return  fittedDen3

end
# sixthOrderFit(n, evsEns4000realizations[:,2], multiPstadev, eta, emax).([1.14 1.33])
# plot(-emax :0.1:emax, sixthOrderFit(n, evsEns4000realizations[:,2], multiPstadev, eta, emax).(collect(-emax :0.1:emax)) )
############################################################################################################################################
# 6th order  unfolding with odd terms, realization by realization
begin
    unfoldedEnergyLevels =  Array{Float64}(UndefInitializer(), size(evsEns4000realizations))
    for k = 1 : nens
        specDenFit = sixthOrderFit(n, evsEns4000realizations[:,k], multiPstadev, eta, emax)
        unfoldedEnergyLevels[:,k] =unfoldTrapz(specDenFit,evsEns4000realizations[:,k])
        if mod(k, 200) == 0
            println(k)
        end
    end 
    open("n=24q=4unfoldedLevels(6thOrderLocalUnfoldingWithOddTerms)0p9Bulk","w") do  ufhquadfile
    writedlm(ufhquadfile,unfoldedEnergyLevels[1:2048*4000])
    end
end

