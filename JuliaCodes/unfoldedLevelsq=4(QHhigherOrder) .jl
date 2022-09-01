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
#specDenq2(x) = 2^(n/2-1)*exp(-x^2/(2*multiPstadev^2))/sqrt(2 * π * multiPstadev^2) #normalize the density so that the integral = total number of states
#specDenq4(x) = 2^(n/2-1) * qHdensityNormalized(x,eta, multiPstadev , 30)  #qHdensityNormalized(x,eta, sig, ntrunc) 
begin
    model0(x, p) = p[2] * qHdensityNormalized.(x,eta, multiPstadev , 30) .* .+(1,  p[1]*qHpolynomial.(x, 6, eta, multiPstadev) )
    p0model0 = [0, 2^(n/2-1)*nens ]
    nbins = 107 #trial and erro so that chi square is minimized
    (xdata, ydata) = histogramData(evsq4Ens4000[:,1], -emax, emax, nbins)

    histoCurveQHLsqFitModel0 = curve_fit(model0, xdata[3: end-2], ydata[3: end-2], p0model0) #exclude edges of the histogram when fitting
    fittedCoefModel0= coef(histoCurveQHLsqFitModel0)
    chiSquareModel0 = sum( histoCurveQHLsqFitModel0.resid.^2 ) #chi square if the density were normalized to 1
    println(chiSquareModel0)
    fittedDenModel0(x) = 2^(n/2-1)* qHdensityNormalized(x,eta, multiPstadev , 30) * (1 + fittedCoefModel0[1]*qHpolynomial(x, 6, eta, multiPstadev) )
end

begin
    model1(x, p) = p[4] * qHdensityNormalized.(x,eta, multiPstadev , 30) .* .+(1,  p[1]*qHpolynomial.(x, 2, eta, multiPstadev),  p[2]*qHpolynomial.(x, 4, eta, multiPstadev) , p[3]*qHpolynomial.(x, 6, eta, multiPstadev) )
    p0model1 = [0, 0, 0, 2^(n/2-1)*nens ]
    nbins = 107
    bulkFraction = 0.99
    (xdata, ydata) = histogramData(evsq4Ens4000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

    histoCurveQHLsqFitModel1 = curve_fit(model1, xdata, ydata, p0model1) #exclude edges of the histogram when fitting
    fittedCoefModel1= coef(histoCurveQHLsqFitModel1)
    chiSquareModel1 = sum( histoCurveQHLsqFitModel1.resid.^2 ) #chi square if the density were normalized to 1
    println(chiSquareModel1)
    fittedDenModel1(x) = 2^(n/2-1)* qHdensityNormalized(x,eta, multiPstadev , 30) * (1 + fittedCoefModel1[1]*qHpolynomial(x, 2, eta, multiPstadev)+ fittedCoefModel1[2]*qHpolynomial(x, 4, eta, multiPstadev) + fittedCoefModel1[3]*qHpolynomial(x, 6, eta, multiPstadev) )
end

begin
    model2(x, p) = p[5] * qHdensityNormalized.(x,eta, multiPstadev , 30) .* .+(1,  p[1]*qHpolynomial.(x, 2, eta, multiPstadev),  p[2]*qHpolynomial.(x, 4, eta, multiPstadev) , p[3]*qHpolynomial.(x, 6, eta, multiPstadev), p[4]*qHpolynomial.(x, 8, eta, multiPstadev) )
    p0model2 = [0, 0, 0, 0, 2^(n/2-1)*nens ]
    nbins = 107
    bulkFraction = 0.99
    (xdata, ydata) = histogramData(evsq4Ens4000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

    histoCurveQHLsqFitModel2 = curve_fit(model2, xdata, ydata, p0model2) #exclude edges of the histogram when fitting
    fittedCoefModel2= coef(histoCurveQHLsqFitModel2)
    chiSquareModel2 = sum( histoCurveQHLsqFitModel2.resid.^2 ) #chi square if the density were normalized to 1
    println(chiSquareModel2)
    fittedDenModel2(x) = 2^(n/2-1)* qHdensityNormalized(x,eta, multiPstadev , 30) * (1 + fittedCoefModel2[1]*qHpolynomial(x, 2, eta, multiPstadev)+ fittedCoefModel2[2]*qHpolynomial(x, 4, eta, multiPstadev) + fittedCoefModel2[3]*qHpolynomial(x, 6, eta, multiPstadev) + fittedCoefModel2[4]*qHpolynomial(x, 8, eta, multiPstadev))
end

begin
    model3(x, p) = p[7] * qHdensityNormalized.(x,eta, multiPstadev , 30) .* .+(1,  p[1]*qHpolynomial.(x, 1, eta, multiPstadev),  p[2]*qHpolynomial.(x, 2, eta, multiPstadev) , p[3]*qHpolynomial.(x, 3, eta, multiPstadev), p[4]*qHpolynomial.(x, 4, eta, multiPstadev) , p[5]*qHpolynomial.(x, 5, eta, multiPstadev),  p[6]*qHpolynomial.(x, 6, eta, multiPstadev) )
    p0model3 = [0, 0, 0, 0,0, 0, 2^(n/2-1)*nens ]
    nbins = 71
    bulkFraction =0.9
    (xdata, ydata) = histogramData(evsq4Ens4000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

    histoCurveQHLsqFitModel3 = curve_fit(model3, xdata, ydata, p0model3) #exclude edges of the histogram when fitting
    fittedCoefModel3= coef(histoCurveQHLsqFitModel3)
    chiSquareModel3 = sum( histoCurveQHLsqFitModel3.resid.^2 ) #chi square if the density were normalized to 1
    println(chiSquareModel3)
    fittedDenModel3(x) = 2^(n/2-1) * qHdensityNormalized.(x,eta, multiPstadev , 30) .* .+(1,  fittedCoefModel3[1]*qHpolynomial.(x, 1, eta, multiPstadev), fittedCoefModel3[2]*qHpolynomial.(x, 2, eta, multiPstadev) , fittedCoefModel3[3]*qHpolynomial.(x, 3, eta, multiPstadev), fittedCoefModel3[4]*qHpolynomial.(x, 4, eta, multiPstadev) , fittedCoefModel3[5]*qHpolynomial.(x, 5, eta, multiPstadev),  fittedCoefModel3[6]*qHpolynomial.(x, 6, eta, multiPstadev) )

end

hquadrature(fittedDenModel3, -emax, emax,reltol=1e-7, abstol=0, maxevals=0)[1]
(xdatafull, ydatafull) = histogramData(evsq4Ens4000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

p1 = plot(xdatafull[1:10], ydatafull[1:10])
fittedOriginalModel3(x)= fittedDenModel3(x)/2^(n/2-1) * fittedCoefModel3[7]
plot!(p1,xdatafull[1:10],fittedOriginalModel3 )
plot(xdatafull[1:10], [ ydatafull[1:10] fittedOriginalModel3.(xdatafull[1:10])], label = ["actual (4000 realizations)" "6th order fitted"], title ="Average spectral density (N=24, q=4)", legend = :bottomright)
savefig("SpectralDenSixthOrderFitEdgeUnderShoot.pdf")
plot(xdatafull, [ ydatafull fittedOriginalModel3.(xdatafull)], label = ["actual (4000 realizations)" "6th order fitted"], title ="Average spectral density (N=24, q=4)", legend = :bottomright)

# begin
#     model4(x, p) = p[8] * qHdensityNormalized.(x,p[7], multiPstadev , 30) .* .+(1,  p[1]*qHpolynomial.(x, 1, p[7],multiPstadev),  p[2]*qHpolynomial.(x, 2, p[7], multiPstadev) , p[3]*qHpolynomial.(x, 3, p[7], multiPstadev), p[4]*qHpolynomial.(x, 4, p[7], multiPstadev) , p[5]*qHpolynomial.(x, 5, p[7], multiPstadev),  p[6]*qHpolynomial.(x, 6, p[7],multiPstadev) )
#     p0model4 = [0, 0, 0,0, 0,0,0, 2^(n/2-1)*nens ]
#     nbins = 105
#     bulkFraction = 0.98
#     (xdata, ydata) = histogramData(evsq4Ens4000[:,1], -emax*bulkFraction, emax*bulkFraction, nbins)

#     histoCurveQHLsqFitModel4 = curve_fit(model4, xdata, ydata, p0model4) #exclude edges of the histogram when fitting
#     fittedCoefModel4= coef(histoCurveQHLsqFitModel4)
#     chiSquareModel4 = sum( histoCurveQHLsqFitModel4.resid.^2 ) #chi square if the density were normalized to 1
#     println(chiSquareModel4)
#     fittedDenModel4(x) = 2^(n/2-1) * qHdensityNormalized.(x,fittedCoefModel4[7], multiPstadev , 30) .* .+(1,  fittedCoefModel4[1]*qHpolynomial.(x, 1, fittedCoefModel4[7], multiPstadev), fittedCoefModel4[2]*qHpolynomial.(x, 2, fittedCoefModel4[7], multiPstadev) , fittedCoefModel4[3]*qHpolynomial.(x, 3, fittedCoefModel4[7], multiPstadev), fittedCoefModel4[4]*qHpolynomial.(x, 4, fittedCoefModel4[7], multiPstadev) , fittedCoefModel4[5]*qHpolynomial.(x, 5, fittedCoefModel4[7], multiPstadev),  fittedCoefModel4[6]*qHpolynomial.(x, 6, fittedCoefModel4[7], multiPstadev) )

# end
# (xdatafull, ydatafull) = histogramData(evsq4Ens4000[:,1], -emax, emax, nbins)
# p1 = plot(xdatafull[1:10], ydatafull[1:10])
# fittedOriginalModel4(x)= fittedDenModel4(x)/2^(n/2-1) * fittedCoefModel4[8]
# plot!(p1,xdatafull[1:10],fittedOriginalModel4 )
# plot(xdatafull[1:10], [ ydatafull[1:10] fittedOriginalModel4.(xdatafull[1:10])], label = ["actual (4000 realizations)" "6th order fitted"], title ="Average spectral density (N=24, q=4)", legend = :bottomright)
# savefig("SpectralDenSixthOrderFitEdgeUnderShoot.pdf")
# plot(xdatafull, [ ydatafull fittedOriginalModel4.(xdatafull)], label = ["actual (4000 realizations)" "6th order fitted"], title ="Average spectral density (N=24, q=4)", legend = :bottomright)

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
############################################################################################################################################
# 6th order ensemble unfolding
unfoldedEnergyLevels =  Array{Float64}(UndefInitializer(), size(evsEns4000realizations))
for k = 1 : nens
    unfoldedEnergyLevels[:,k] =unfoldTrapz(fittedDenModel1,evsEns4000realizations[:,k])
    if mod(k, 200) == 0
        println(k)
    end
end 
open("n=24q=4unfoldedLevels(6thOrderEnsembleUnfolding)0p99Bulk","w") do  ufhquadfile
writedlm(ufhquadfile,unfoldedEnergyLevels[1:2048*4000])
end
############################################################################################################################################
# 8th order ensemble unfolding
unfoldedEnergyLevels =  Array{Float64}(UndefInitializer(), size(evsEns4000realizations))
for k = 1 : nens
    unfoldedEnergyLevels[:,k] =unfoldTrapz(fittedDenModel2,evsEns4000realizations[:,k])
    if mod(k, 200) == 0
        println(k)
    end
end 
open("n=24q=4unfoldedLevels(8thOrderEnsembleUnfolding)0p99bulk","w") do  ufhquadfile
writedlm(ufhquadfile,unfoldedEnergyLevels[1:2048*4000])
end
############################################################################################################################################
# 6th order ensemble unfolding with odd terms
begin
    unfoldedEnergyLevels =  Array{Float64}(UndefInitializer(), size(evsEns4000realizations))
    for k = 1 : nens
        unfoldedEnergyLevels[:,k] =unfoldTrapz(fittedDenModel3,evsEns4000realizations[:,k])
        if mod(k, 200) == 0
            println(k)
        end
    end 
    open("n=24q=4unfoldedLevels(6thOrderEnsembleUnfoldingWithOddTerms)0p9Bulk","w") do  ufhquadfile
    writedlm(ufhquadfile,unfoldedEnergyLevels[1:2048*4000])
    end
end

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


