
function  sffConnectedCompute(times::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, enRealizations::Array{Float64,2})
    #compute the connected spectral form factor from given realizations of spectra (unfolded or not)
    #times: time range and step length, e.g. times = 0 : 0.1 : 5
    #enRealizations: realizations of spectra, each column is one realizations
    sffconnected = zeros(length(times))

    for k = 1: length(times)
    allphases = exp.(im * times[k] * enRealizations)
    phaseColumnSum= sum(eachrow(allphases))
    sffconnected[k] =  norm(phaseColumnSum)^2/nens -abs(sum(phaseColumnSum))^2/nens^2
    end

    return [collect(times) sffconnected/length(enRealizations[:,1])] #normalize late-time sff to 1
end