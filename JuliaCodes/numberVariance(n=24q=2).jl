begin
    include("/home/yiyang/Desktop/Codes/julia/SYK/numberVarianceCompute.jl")
    include("/home/yiyang/Desktop/Codes/julia/SYK/RMTfunctions.jl")
    include("/home/yiyang/Desktop/Codes/julia/SYK/unfoldingMethods.jl")

    using LinearAlgebra
    using DelimitedFiles
    using Plots
    import StatsPlots 
    import PyPlot
    using LaTeXStrings
    using StatsBase
    using LsqFit
end

n = 24
q =2 
nens = 4000


evsq2Ens4000 = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=2unfoldedLevels(hquadMethod)"), (2^(Int(n/2)-1), nens) )
intervalMaxLenghFraction = 1/2
intervalMaxLength = Int( length(evsq2Ens4000[:, 1]) * intervalMaxLenghFraction )
samplingPoints =  0.01 : 0.1 : intervalMaxLength
numberVariance = zeros( length(samplingPoints) )

for k = 1 : length(samplingPoints)
    intervalLength = samplingPoints[k]
    for j = 1: nens
        numberVariance[k] = numberVariance[k] + nVcomputeNoOverlap(evsq2Ens4000[:, j], intervalLength)
    end
    numberVariance[k] = numberVariance[k]/nens
end

open("numberVariancen=24q=2NoOverlap","w") do file
    writedlm(file, [collect(samplingPoints) numberVariance] )
end


samplingPoints = readdlm("/home/yiyang/Desktop/Codes/julia/numberVariancen=24q=2NoOverlap")[:,1]
numberVariance =  readdlm("/home/yiyang/Desktop/Codes/julia/numberVariancen=24q=2NoOverlap")[:,2]
plot(collect(samplingPoints), numberVariance, label= nothing)
plot(collect(samplingPoints)[1:400], numberVariance[1:400], label= nothing)
######################################################################################################################################################################################################################################
# sixth order ensemble unfolding

n = 24
q =2 
nens = 8000
evsq2Ens8000sixthOrder = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=2unfoldedLevelsEns8000(6thOrderEnsembleUnfoldingWithOddTerms)0p9Bulk"), (2^(Int(n/2)-1), nens) )
# StatsPlots.histogram(evsq4Ens4000sixthOrder[1: 2048*4000])
# (xdata, ydata)= histogramData(evsq4Ens4000sixthOrder[1: 2048*4000], 0, 2048, 100)
# plot(xdata, ydata)
begin
    xsampling = 0 : 1 : length(evsq2Ens8000sixthOrder[:, 1])
    stairs = zeros(length(xsampling))
    for k = 1 :  length(xsampling)
        x = xsampling[k] 
        for j = 1:nens
            stairs[k] = stairs[k] + stairCase(x, evsq2Ens8000sixthOrder[:, j]) - x
        end
        stairs[k] = stairs[k]/nens
    end

    pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", title="Sixth-order ensemble unfolding with odd terms", legend =:topright)
end

savefig("unfoldedStaircase(N=24,q=2)AveSixthOrderWithOddTermsFit.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best, title = "")
savefig("unfoldedStaircase(N=24,q=2)Realizations.pdf")

begin
    intervalMaxLenghFraction = 1/4
    intervalMaxLength = round(Int, length(evsq2Ens8000sixthOrder[:, 1]) * intervalMaxLenghFraction )
    samplingPoints =  0.01 : 0.1 : intervalMaxLength
    numberVariance = zeros( length(samplingPoints) )
    for k = 1 : length(samplingPoints)
        intervalLength = samplingPoints[k]
        for j = 1: nens
            numberVariance[k] = numberVariance[k] + nVcomputeNoOverlap(evsq2Ens8000sixthOrder[:, j], intervalLength)
        end
        numberVariance[k] = numberVariance[k]/nens
    end

    plotNV = plot(samplingPoints, numberVariance)
end

open("numberVarianceN24q2SixthOrderEnsembleUnfolding","w") do ensUnf
    writedlm(ensUnf, [collect(samplingPoints)  numberVariance ])
end
plotNV50 = plot(samplingPoints[1:500], [numberVariance[1:500] nVGOE.(samplingPoints[1:500]) samplingPoints[1:500]], label =["N=24, q=2 SYK, 6th order ensemble unfolding" "GOE" "Poisson"], legend=:topleft)
plotNV50 = plot(samplingPoints, numberVariance , label =["N=24, q=2 SYK, 6th order ensemble unfolding" ], legend=:topleft)
plotNV50 = plot(samplingPoints[1:20], [numberVariance[1:20] nVGOE.(samplingPoints[1:20])], label =["N=24, q=2 SYK, 6th order ensemble unfolding" "GOE"], legend=:topleft)
savefig("numberVariance6thOrderEnsUnfolding(N=24,q=2).pdf")
plotNV50 = plot(samplingPoints[1:4000], [numberVariance[1:4000] samplingPoints[1:4000]], label =["N=24, q=2 SYK, 6th order ensemble unfolding - Poisson"], legend=:topleft)

modelNV(x, p)=p[1].+ p[2] * x  .+ p[3]* x.^2   
p0modelNV = [1, 0 , 1/2/binomial(n,q)]
longRangeNVFit = curve_fit(modelNV, samplingPoints[150:end], numberVariance[150: end]-  nVGOE.(samplingPoints[150:end]) , p0modelNV)
longRangeCoe = coef(longRangeNVFit)

longRangeFunction(x) = longRangeCoe[1] + longRangeCoe[2]*x + longRangeCoe[3]*x^2
plotLongRange =  plot(samplingPoints, [numberVariance   nVGOE.(samplingPoints) nVGOE.(samplingPoints)+longRangeFunction.(samplingPoints)], linewidth = [1 1.5 1.5], title="N=24, q=4 SYK number variance", xlabel = L"\bar n", ylabel=L"\Sigma^2(\bar n)",label =["6th order ensemble unfolding" "GOE" L"GOE + 0.126 -0.002\bar n + 4.634\times 10^{-5}\bar n^2 "], legend=:topleft)
savefig("numberVarianceLongRange(N=24,q=4).pdf")

# longRangeNVFitRMTsubtracted = curve_fit(modelNV, samplingPoints[1:5000], numberVariance[1: 5000] - nVGOE.(samplingPoints[1:5000]), p0modelNV)
# longRangeCoeRMTsubtracted = coef(longRangeNVFitRMTsubtracted)
# plot(samplingPoints[1:6000], numberVariance[1: 6000] - nVGOE.(samplingPoints[1:6000]))

####################################################################################################################################################################################################################################
# sixth order local unfolding

n = 24
q =2 
nens = 8000
evsq2Ens8000sixthOrder = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=2unfoldedLevelsEns8000(6thOrderEnsembleUnfoldingWithOddTerms)0p9Bulk"), (2^(Int(n/2)-1), nens) )
# StatsPlots.histogram(evsq4Ens4000sixthOrder[1: 2048*4000])
# (xdata, ydata)= histogramData(evsq4Ens4000sixthOrder[1: 2048*4000], 0, 2048, 100)
# plot(xdata, ydata)
begin
    xsampling = 0 : 1 : length(evsq2Ens8000sixthOrder[:, 1])
    stairs = zeros(length(xsampling))
    for k = 1 :  length(xsampling)
        x = xsampling[k] 
        for j = 1:nens
            stairs[k] = stairs[k] + stairCase(x, evsq2Ens8000sixthOrder[:, j]) - x
        end
        stairs[k] = stairs[k]/nens
    end

    pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", title="Sixth-order ensemble unfolding with odd terms", legend =:topright)
end

savefig("unfoldedStaircase(N=24,q=2)AveSixthOrderWithOddTermsFit.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best, title = "")
savefig("unfoldedStaircase(N=24,q=2)Realizations.pdf")

begin
    intervalMaxLenghFraction = 1/4
    intervalMaxLength = round(Int, length(evsq2Ens8000sixthOrder[:, 1]) * intervalMaxLenghFraction )
    samplingPoints =  0.01 : 0.1 : intervalMaxLength
    numberVariance = zeros( length(samplingPoints) )
    for k = 1 : length(samplingPoints)
        intervalLength = samplingPoints[k]
        for j = 1: nens
            numberVariance[k] = numberVariance[k] + nVcomputeNoOverlap(evsq2Ens8000sixthOrder[:, j], intervalLength)
        end
        numberVariance[k] = numberVariance[k]/nens
    end

    plotNV = plot(samplingPoints, numberVariance)
end

open("numberVarianceN24q2SixthOrderLocalUnfolding","w") do ensUnf
    writedlm(ensUnf, [collect(samplingPoints)  numberVariance ])
end
plotNV50 = plot(samplingPoints[1:20], [numberVariance[1:20] nVGOE.(samplingPoints[1:20]) samplingPoints[1:20]], label =["N=24, q=2 SYK, 6th order local unfolding" "GOE" "Poisson"], legend=:topleft)
plotNV50 = plot(samplingPoints, numberVariance , label =["N=24, q=2 SYK, 6th order local unfolding" ], legend=:topleft)
plotNV50 = plot(samplingPoints[1:20], [numberVariance[1:20] nVGOE.(samplingPoints[1:20])], label =["N=24, q=2 SYK, 6th order ensemble unfolding" "GOE"], legend=:topleft)
savefig("numberVariance6thOrderLocalUnfolding(N=24,q=2)ZoomIn.pdf")
plotNV50 = plot(samplingPoints[1:4000], [numberVariance[1:4000] samplingPoints[1:4000]], label =["N=24, q=2 SYK, 6th order ensemble unfolding - Poisson"], legend=:topleft)

modelNV(x, p)=p[1].+ p[2] * x  .+ p[3]* x.^2   
p0modelNV = [1, 0 , 1/2/binomial(n,q)]
longRangeNVFit = curve_fit(modelNV, samplingPoints[150:end], numberVariance[150: end]-  nVGOE.(samplingPoints[150:end]) , p0modelNV)
longRangeCoe = coef(longRangeNVFit)

longRangeFunction(x) = longRangeCoe[1] + longRangeCoe[2]*x + longRangeCoe[3]*x^2
plotLongRange =  plot(samplingPoints, [numberVariance   nVGOE.(samplingPoints) nVGOE.(samplingPoints)+longRangeFunction.(samplingPoints)], linewidth = [1 1.5 1.5], title="N=24, q=4 SYK number variance", xlabel = L"\bar n", ylabel=L"\Sigma^2(\bar n)",label =["6th order ensemble unfolding" "GOE" L"GOE + 0.126 -0.002\bar n + 4.634\times 10^{-5}\bar n^2 "], legend=:topleft)
savefig("numberVarianceLongRange(N=24,q=4).pdf")

# longRangeNVFitRMTsubtracted = curve_fit(modelNV, samplingPoints[1:5000], numberVariance[1: 5000] - nVGOE.(samplingPoints[1:5000]), p0modelNV)
# longRangeCoeRMTsubtracted = coef(longRangeNVFitRMTsubtracted)
# plot(samplingPoints[1:6000], numberVariance[1: 6000] - nVGOE.(samplingPoints[1:6000]))

####################################################################################################################################################################################################################################