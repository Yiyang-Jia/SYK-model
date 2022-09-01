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
q =4
nens = 4000


evsq4Ens4000 = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=4unfoldedLevels(trapezoidMethod)"), (2^(Int(n/2)-1), nens) )
intervalMaxLenghFraction = 1/2
intervalMaxLength = Int( length(evsq4Ens4000[:, 1]) * intervalMaxLenghFraction )
samplingPoints =  0.01 : 0.1 : intervalMaxLength
numberVariance = zeros( length(samplingPoints) )
########################################################################################################################################################################
for k = 1 : length(samplingPoints)
    intervalLength = samplingPoints[k]
    overalapOnsetThreshold = 1/100 #length(evsq4Ens4000[:, 1]) * overalapOnsetThreshold is the interval length after which we start to introduce interval overalaps
    if intervalLength < length(evsq4Ens4000[:, 1]) * overalapOnsetThreshold
        overlapFrac = 0.0
        for j = 1: nens
        numberVariance[k] = numberVariance[k] + nVcompute(evsq4Ens4000[:, j], intervalLength, overlapFrac)
        end
        numberVariance[k] = numberVariance[k]/nens
    elseif  length(evsq4Ens4000[:, 1]) * overalapOnsetThreshold < intervalLength < length(evsq4Ens4000[:, 1]) * 0.33
        # overlapFrac =   0.5 * ( intervalLength/length(evsq4Ens4000[:, 1]) -overalapOnsetThreshold ) / (intervalMaxLenghFraction - overalapOnsetThreshold) 
        overlapFrac =0.2
        for j = 1: nens
            numberVariance[k] = numberVariance[k] + nVcompute(evsq4Ens4000[:, j], intervalLength, overlapFrac)
        end
        numberVariance[k] = numberVariance[k]/nens
    else 
        # overlapFrac =   0.5 * ( intervalLength/length(evsq4Ens4000[:, 1]) -overalapOnsetThreshold ) / (intervalMaxLenghFraction - overalapOnsetThreshold) 
        overlapFrac =0.4
        for j = 1: nens
            numberVariance[k] = numberVariance[k] + nVcompute(evsq4Ens4000[:, j], intervalLength, overlapFrac)
        end
        numberVariance[k] = numberVariance[k]/nens
    end
end

open("numberVariancen=24q=4StepOverlapThresholds(Problematic)","w") do file
    writedlm(file, [collect(samplingPoints) numberVariance] )
end
plot(collect(samplingPoints), numberVariance, label= nothing)
plot(collect(samplingPoints)[1:400], numberVariance[1:400], label= nothing)
###################################################################################################################################################
for k = 1 : length(samplingPoints)
    intervalLength = samplingPoints[k]
    for j = 1: nens
        numberVariance[k] = numberVariance[k] + nVcomputeNoOverlap(evsq4Ens4000[:, j], intervalLength)
    end
    numberVariance[k] = numberVariance[k]/nens
end

plotNV = plot(samplingPoints, numberVariance)
quad(m) = 1/2 * 1/binomial(n, q) * m^2
plot!(plotNV, 0:0.1:1000, quad)
plotNVsmall = plot(samplingPoints[1:500], numberVariance[1:500])
plot!(plotNVsmall, 0:0.1:50, nVGOE)
open("numberVariancen=24q=4NoOverlap","w") do file
    writedlm(file, [collect(samplingPoints) numberVariance] )
end


# samplingPointsq2 = readdlm("/home/yiyang/Desktop/Codes/julia/numberVariancen=24q=2NoOverlap")[:,1]
# numberVarianceq2 =  readdlm("/home/yiyang/Desktop/Codes/julia/numberVariancen=24q=2NoOverlap")[:,2]

# samplingPointsq4 = readdlm("/home/yiyang/Desktop/Codes/julia/numberVariancen=24q=4NoOverlap")[:,1]
# numberVarianceq4 =  readdlm("/home/yiyang/Desktop/Codes/julia/numberVariancen=24q=4NoOverlap")[:,2]
# plot(collect(samplingPointsq4), numberVarianceq4, label= nothing)
# nvgoe = plot(0:0.1:40, nVGOE, label = "GOE", title = "Number variance", xlabel = L"$\bar{n}$", ylabel = L"\Sigma^2(\bar n)")
# plot!(nvgoe, collect(samplingPointsq4)[1:400],   numberVarianceq4[1:400], label= "N=24, q=4 SYK")
# savefig("numberVariance(N=24q=4).pdf")

# quad(m) = 1/2 * 1/binomial(n, q) * m^2
# quadplot = plot(0:0.1:300, quad, title = "Number variance (N=24, q=4)", xlabel = L"$\bar{n}$", ylabel = L"\Sigma^2(\bar n)", label = L"{\bar{n}^2}/{2\binom{N}{q}}")
# plot!(quadplot, collect(samplingPointsq4)[1:3000], numberVarianceq4[1:3000], label= "SYK data")
# savefig("numberVariance(N=24q=4)LargeSeparation.pdf")

# rmtPlot= plot(0.0001:0.01:1, [nVGOE, nVGUE ,nVGSE])
# plot!(rmtPlot, samplingPointsq2[1:10], numberVarianceq2[1:10])





println(nVcompute(evsq4Ens4000[:, 3000], 100, 0))

for k = 1 : length(samplingPoints)
    intervalLength = samplingPoints[k]
    for j = 1: nens
        numberVariance[k] = numberVariance[k] + nVcompute(evsq4Ens4000[:, j], intervalLength, 0)
    end
    numberVariance[k] = numberVariance[k]/nens
end
plot(samplingPoints, numberVariance, xlabel = "interval of length n", ylabel = L"\bar n", label = nothing)
savefig("NumberVSintervalLength(N=24,q=4).pdf")



xsampling = 0 : 1 : length(evsq4Ens4000[:, 1])
stairs = zeros(length(xsampling))
for k = 1 :  length(xsampling)
    x = xsampling[k] 
    for j = 1:nens
        stairs[k] = stairs[k] + stairCase(x, evsq4Ens4000[:, j]) - x
    end
    stairs[k] = stairs[k]/nens
end

pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", legend =:bottomright)
savefig("unfoldedStaircase(N=24,q=4)Ave.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best)
savefig("unfoldedStaircase(N=24,q=4)Realizations.pdf")
####################################################################################################################################################################################################################################
# 6th order QH ensemble unfolding
evsq4Ens4000sixthOrder = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=2unfoldedLevels(6thOrderEnsembleUnfolding)"), (2^(Int(n/2)-1), nens) )

xsampling = 0 : 1 : length(evsq4Ens4000sixthOrder[:, 1])
stairs = zeros(length(xsampling))
for k = 1 :  length(xsampling)
    x = xsampling[k] 
    for j = 1:nens
        stairs[k] = stairs[k] + stairCase(x, evsq4Ens4000sixthOrder[:, j]) - x
    end
    stairs[k] = stairs[k]/nens
end

pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", title="Sixth-order ensemble unfolding", legend =:topright)
savefig("unfoldedStaircase(N=24,q=4)AveSixthOrderFit.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best)
savefig("unfoldedStaircase(N=24,q=4)Realizations.pdf")

####################################################################################################################################################################################################################################
# 8th order QH ensemble unfolding
evsq4Ens4000eighthOrder = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=2unfoldedLevels(8thOrderEnsembleUnfolding)"), (2^(Int(n/2)-1), nens) )

xsampling = 0 : 1 : length(evsq4Ens4000eighthOrder[:, 1])
stairs = zeros(length(xsampling))
for k = 1 :  length(xsampling)
    x = xsampling[k] 
    for j = 1:nens
        stairs[k] = stairs[k] + stairCase(x, evsq4Ens4000eighthOrder[:, j]) - x
    end
    stairs[k] = stairs[k]/nens
end

pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", title="Eighth-order ensemble unfolding", legend =:topright)
savefig("unfoldedStaircase(N=24,q=4)AveSixthOrderFit.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best)
savefig("unfoldedStaircase(N=24,q=4)Realizations.pdf")

# another 8th order QH ensemble unfolding
evsq4Ens4000eighthOrder = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=2unfoldedLevels(8thOrderEnsembleUnfolding)0p99bulk"), (2^(Int(n/2)-1), nens) )

xsampling = 0 : 1 : length(evsq4Ens4000eighthOrder[:, 1])
stairs = zeros(length(xsampling))
for k = 1 :  length(xsampling)
    x = xsampling[k] 
    for j = 1:nens
        stairs[k] = stairs[k] + stairCase(x, evsq4Ens4000eighthOrder[:, j]) - x
    end
    stairs[k] = stairs[k]/nens
end

pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", title="Eighth-order ensemble unfolding 99% bulk", legend =:topright)
savefig("unfoldedStaircase(N=24,q=4)AveSixthOrderFit.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best)
savefig("unfoldedStaircase(N=24,q=4)Realizations.pdf")

####################################################################################################################################################################################################################################
# 6th order QH ensemble unfolding with odd terms, 0.98*emax
evsq4Ens4000sixthOrder = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=2unfoldedLevels(6thOrderEnsembleUnfoldingWithOddTerms)0p98Bulk"), (2^(Int(n/2)-1), nens) )

xsampling = 1 : 1 : length(evsq4Ens4000sixthOrder[:, 1])
stairs = zeros(length(xsampling))
for k = 1 :  length(xsampling)
    x = xsampling[k] 
    for j = 1:nens
        stairs[k] = stairs[k] + stairCase(x, evsq4Ens4000sixthOrder[:, j]) - x
    end
    stairs[k] = stairs[k]/nens
end

pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", title="Sixth-order ensemble unfolding with odd terms", legend =:topright)
savefig("unfoldedStaircase(N=24,q=4)AveSixthOrderWithOddTermsFit.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best, title = "")
savefig("unfoldedStaircase(N=24,q=4)Realizations.pdf")


####################################################################################################################################################################################################################################
# 6th order QH ensemble unfolding with odd terms, 0.9*emax
evsq4Ens4000sixthOrder = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=4unfoldedLevels(6thOrderEnsembleUnfoldingWithOddTerms)0p9Bulk"), (2^(Int(n/2)-1), nens) )
# StatsPlots.histogram(evsq4Ens4000sixthOrder[1: 2048*4000])
# (xdata, ydata)= histogramData(evsq4Ens4000sixthOrder[1: 2048*4000], 0, 2048, 100)
# plot(xdata, ydata)
begin
    xsampling = 0 : 1 : length(evsq4Ens4000sixthOrder[:, 1])
    stairs = zeros(length(xsampling))
    for k = 1 :  length(xsampling)
        x = xsampling[k] 
        for j = 1:nens
            stairs[k] = stairs[k] + stairCase(x, evsq4Ens4000sixthOrder[:, j]) - x
        end
        stairs[k] = stairs[k]/nens
    end

    pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", title="Sixth-order ensemble unfolding with odd terms", legend =:topright)
end

savefig("unfoldedStaircase(N=24,q=4)AveSixthOrderWithOddTermsFit.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best, title = "")
savefig("unfoldedStaircase(N=24,q=4)Realizations.pdf")

begin
    intervalMaxLenghFraction = 1/4
    intervalMaxLength = round(Int, length(evsq4Ens4000sixthOrder[:, 1]) * intervalMaxLenghFraction )
    samplingPoints =  0.01 : 0.1 : intervalMaxLength
    numberVariance = zeros( length(samplingPoints) )
    for k = 1 : length(samplingPoints)
        intervalLength = samplingPoints[k]
        for j = 1: nens
            numberVariance[k] = numberVariance[k] + nVcomputeNoOverlap(evsq4Ens4000sixthOrder[:, j], intervalLength)
        end
        numberVariance[k] = numberVariance[k]/nens
    end

    plotNV = plot(samplingPoints, numberVariance)
end

open("numberVarianceN24q4SixthOrderEnsembleUnfolding","w") do ensUnf
    writedlm(ensUnf, [collect(samplingPoints)  numberVariance ])
end
plotNV50 = plot(samplingPoints[1:500], [numberVariance[1:500] nVGOE.(samplingPoints[1:500])], label =["N=24, q=4 SYK, 6th order ensemble unfolding" "GOE"], legend=:bottomright)

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
# 6th order QH local unfolding with odd terms, 0.9*emax
evsq4Ens4000sixthOrder = reshape(readdlm("/home/yiyang/Desktop/Codes/julia/n=24q=4unfoldedLevels(6thOrderLocalUnfoldingWithOddTerms)0p9Bulk"), (2^(Int(n/2)-1), nens) )
StatsPlots.histogram(evsq4Ens4000sixthOrder[1: 2048*4000])
(xdata, ydata)= histogramData(evsq4Ens4000sixthOrder[1: 2048*4000], 0, 2048, 1000)
plot(xdata, ydata)
begin
    xsampling = 0 : 1 : length(evsq4Ens4000sixthOrder[:, 1])
    stairs = zeros(length(xsampling))
    for k = 1 :  length(xsampling)
        x = xsampling[k] 
        for j = 1:nens
            stairs[k] = stairs[k] + stairCase(x, evsq4Ens4000sixthOrder[:, j]) - x
        end
        stairs[k] = stairs[k]/nens
    end

    pAve = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E)- E", label = "Ensemble average", title="Sixth-order ensemble unfolding with odd terms", legend =:topright)
end

savefig("unfoldedStaircase(N=24,q=4)AveSixthOrderWithOddTermsFit.pdf")
p2 = plot(xsampling, stairs, xlabel = L"E  \ (unfolded)", ylabel = L" n(E) - E", label = "2nd realization")
p11 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "11th realization")
p300 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "300th realization")
p3099 = plot(xsampling, stairs, xlabel = L"E \ (unfolded)", ylabel = L" n(E) - E", label = "3099th realization", legend= :best)
plot(p2,p11,p300,p3099, layout=4, legend =:best, title = "")
savefig("unfoldedStaircase(N=24,q=4)Realizations.pdf")

begin
    intervalMaxLenghFraction = 1/3
    intervalMaxLength = round(Int, length(evsq4Ens4000sixthOrder[:, 1]) * intervalMaxLenghFraction )
    samplingPoints =  0.01 : 0.1 : intervalMaxLength
    numberVariance = zeros( length(samplingPoints) )
    for k = 1 : length(samplingPoints)
        intervalLength = samplingPoints[k]
        for j = 1: nens
            numberVariance[k] = numberVariance[k] + nVcomputeNoOverlap(evsq4Ens4000sixthOrder[:, j], intervalLength)
        end
        numberVariance[k] = numberVariance[k]/nens
    end

    plotNV = plot(samplingPoints, numberVariance)
end
open("numberVarianceN24q4SixthOrderLocalUnfolding","w") do localUnf
    writedlm(localUnf, [collect(samplingPoints)  numberVariance ])
end


plotNV200 = plot(samplingPoints[1:2000], [numberVariance[1:2000] nVGOE.(samplingPoints[1:2000])], label =["N=24, q=4 SYK, 6th order local unfolding" "GOE"], legend=:bottomright)
plotNV600 = plot(samplingPoints[1:6000], [numberVariance[1:6000] nVGOE.(samplingPoints[1:6000])], label =["N=24, q=4 SYK, 6th order local unfolding" "GOE"], legend=:bottomright)

# modelNV(x, p)=p[1].+ p[2] * x  .+ p[3]* x.^2   
# p0modelNV = [1, 0 , 1/2/binomial(n,q)]
# longRangeNVFit = curve_fit(modelNV, samplingPoints[100:6800], numberVariance[100: 6800] , p0modelNV)
# longRangeCoe = coef(longRangeNVFit)

# longRangeFunction(x) = longRangeCoe[1] + longRangeCoe[2]*x + longRangeCoe[3]*x^2
# plotLongRange =  plot(samplingPoints, [numberVariance longRangeFunction.(samplingPoints)], linewidth = [1 1.5 ], title="N=24, q=4 SYK number variance", xlabel = L"\bar n", ylabel=L"\Sigma^2(\bar n)",label =["6th order ensemble unfolding" L"1.3823 -0.0015\bar n + 4.74 \times 10^{-5}\bar n^2 "], legend=:bottomright)
savefig("numberVariance(N=24q=4)6thOrderLocalUnfolding.pdf")
