begin
include("/home/yiyang/Desktop/Codes/julia/SYK/unfoldingMethods.jl")
include("/home/yiyang/Desktop/Codes/julia/SYK/qHermite.jl")
using LinearAlgebra
using DelimitedFiles
using Plots
using StatsBase
using LsqFit
import StatsPlots 
import PyPlot
end


n = 24
q = 2
nens = 8000
evsq2Ens4000 = readdlm("/home/yiyang/Desktop/Codes/julia/n=24+q=2+ens=4000.txt")
evsEns4000MATLAB = readdlm("n=24q=2eigenens4000fromMATLAB.txt")
evsEns8000 = [evsq2Ens4000;evsEns4000MATLAB]
evsEns8000realizations =  reshape(evsEns8000,(2^(Int(n/2)-1),nens))
#StatsPlots.histogram(evsq4Ens4000)

stadev = sqrt(factorial(q-1)/(2^q * n^(q-1)))
multiPstadev = sqrt( binomial(n,q) * stadev^2 )
eta = suppressionFactor(n, q) #"SYK/qHermite.jl"
emax = groundStateE(multiPstadev, eta)
########################################################################################################################################################################################################################################################################################
# 6th order ensemble unfolding with odd terms

begin
    model3(x, p) = p[7] * qHdensityNormalized.(x,eta, multiPstadev , 30) .* .+(1,  p[1]*qHpolynomial.(x, 1, eta, multiPstadev),  p[2]*qHpolynomial.(x, 2, eta, multiPstadev) , p[3]*qHpolynomial.(x, 3, eta, multiPstadev), p[4]*qHpolynomial.(x, 4, eta, multiPstadev) , p[5]*qHpolynomial.(x, 5, eta, multiPstadev),  p[6]*qHpolynomial.(x, 6, eta, multiPstadev) )
    p0model3 = [0, 0, 0, 0,0, 0, 2^(n/2-1)*nens ]
    nbins = 110
    bulkFraction =0.9
    (xdata, ydata) = histogramData(evsEns8000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

    histoCurveQHLsqFitModel3 = curve_fit(model3, xdata, ydata, p0model3) #exclude edges of the histogram when fitting
    fittedCoefModel3= coef(histoCurveQHLsqFitModel3)
    chiSquareModel3 = sum( histoCurveQHLsqFitModel3.resid.^2 ) #chi square if the density were normalized to 1
    println(chiSquareModel3)
    fittedDenModel3(x) = 2^(n/2-1) * qHdensityNormalized(x,eta, multiPstadev , 30) * +(1,  fittedCoefModel3[1]*qHpolynomial(x, 1, eta, multiPstadev), fittedCoefModel3[2]*qHpolynomial(x, 2, eta, multiPstadev) , fittedCoefModel3[3]*qHpolynomial(x, 3, eta, multiPstadev), fittedCoefModel3[4]*qHpolynomial(x, 4, eta, multiPstadev) , fittedCoefModel3[5]*qHpolynomial(x, 5, eta, multiPstadev),  fittedCoefModel3[6]*qHpolynomial(x, 6, eta, multiPstadev) )

end

hquadrature(fittedDenModel3, -emax, emax,reltol=1e-7, abstol=0, maxevals=0)[1]
(xdatafull, ydatafull) = histogramData(evsEns8000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

p1 = plot(xdatafull[1:10], ydatafull[1:10])
fittedOriginalModel3(x)= fittedDenModel3(x)/2^(n/2-1) * fittedCoefModel3[7]
plot!(p1,xdatafull[1:10],fittedOriginalModel3 )
plot(xdatafull[1:10], [ ydatafull[1:10] fittedOriginalModel3.(xdatafull[1:10])], label = ["actual (8000 realizations)" "6th order fitted"], title ="Average spectral density (N=24, q=2)", legend = :bottomright)
savefig("SpectralDenSixthOrderFitEdgeUnderShoot.pdf")
plot(xdatafull, [ ydatafull fittedOriginalModel3.(xdatafull)], label = ["actual (4000 realizations)" "6th order fitted"], title ="Average spectral density (N=24, q=2)", legend = :bottomright)


begin
    unfoldedEnergyLevels =  Array{Float64}(UndefInitializer(), size(evsEns8000realizations))
    for k = 1 : nens
        unfoldedEnergyLevels[:,k] =unfoldTrapz(fittedDenModel3,evsEns8000realizations[:,k])
        if mod(k, 200) == 0
            println(k)
        end
    end 
    open("n=24q=2unfoldedLevelsEns8000(6thOrderEnsembleUnfoldingWithOddTerms)0p9Bulk","w") do  ufhquadfile
    writedlm(ufhquadfile,unfoldedEnergyLevels[1:2^(Int(n/2)-1)*nens])
    end
end


############################################################################################################################################################################################################################################################################################################################################################################
# 8th order ensemble unfolding with odd terms

begin
    model4(x, p) = p[9] * qHdensityNormalized.(x,eta, multiPstadev , 30) .* .+(1,  p[1]*qHpolynomial.(x, 1, eta, multiPstadev),  p[2]*qHpolynomial.(x, 2, eta, multiPstadev) , p[3]*qHpolynomial.(x, 3, eta, multiPstadev), p[4]*qHpolynomial.(x, 4, eta, multiPstadev) , p[5]*qHpolynomial.(x, 5, eta, multiPstadev),  p[6]*qHpolynomial.(x, 6, eta, multiPstadev),  p[7]*qHpolynomial.(x, 7, eta, multiPstadev) ,  p[8]*qHpolynomial.(x, 8, eta, multiPstadev)  )
    p0model4 = [0, 0, 0, 0,0, 0, 0,0, 2^(n/2-1)*nens ]
    nbins = 110
    bulkFraction =0.9
    (xdata, ydata) = histogramData(evsEns8000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

    histoCurveQHLsqFitModel4 = curve_fit(model4, xdata, ydata, p0model4) #exclude edges of the histogram when fitting
    fittedCoefModel4= coef(histoCurveQHLsqFitModel4)
    chiSquareModel4 = sum( histoCurveQHLsqFitModel4.resid.^2 ) #chi square if the density were normalized to 1
    println(chiSquareModel4)
    fittedDenModel4(x) = 2^(n/2-1) * qHdensityNormalized(x,eta, multiPstadev , 30) * +(1,  fittedCoefModel4[1]*qHpolynomial(x, 1, eta, multiPstadev), fittedCoefModel4[2]*qHpolynomial(x, 2, eta, multiPstadev) , fittedCoefModel4[3]*qHpolynomial(x, 3, eta, multiPstadev), fittedCoefModel4[4]*qHpolynomial(x, 4, eta, multiPstadev) , fittedCoefModel4[5]*qHpolynomial(x, 5, eta, multiPstadev),  fittedCoefModel4[6]*qHpolynomial(x, 6, eta, multiPstadev),  fittedCoefModel4[7]*qHpolynomial(x, 7, eta, multiPstadev),  fittedCoefModel4[8]*qHpolynomial(x, 8, eta, multiPstadev) )

end

hquadrature(fittedDenModel3, -emax, emax,reltol=1e-7, abstol=0, maxevals=0)[1]
(xdatafull, ydatafull) = histogramData(evsEns8000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

p1 = plot(xdatafull[1:10], ydatafull[1:10])
fittedOriginalModel3(x)= fittedDenModel3(x)/2^(n/2-1) * fittedCoefModel3[7]
plot!(p1,xdatafull[1:10],fittedOriginalModel3 )
plot(xdatafull[1:10], [ ydatafull[1:10] fittedOriginalModel3.(xdatafull[1:10])], label = ["actual (8000 realizations)" "6th order fitted"], title ="Average spectral density (N=24, q=2)", legend = :bottomright)
savefig("SpectralDenSixthOrderFitEdgeUnderShoot.pdf")
plot(xdatafull, [ ydatafull fittedOriginalModel3.(xdatafull)], label = ["actual (4000 realizations)" "6th order fitted"], title ="Average spectral density (N=24, q=2)", legend = :bottomright)


begin
    unfoldedEnergyLevels =  Array{Float64}(UndefInitializer(), size(evsEns8000realizations))
    for k = 1 : nens
        unfoldedEnergyLevels[:,k] =unfoldTrapz(fittedDenModel3,evsEns8000realizations[:,k])
        if mod(k, 200) == 0
            println(k)
        end
    end 
    open("n=24q=2unfoldedLevelsEns8000(6thOrderEnsembleUnfoldingWithOddTerms)0p9Bulk","w") do  ufhquadfile
    writedlm(ufhquadfile,unfoldedEnergyLevels[1:2^(Int(n/2)-1)*nens])
    end
end