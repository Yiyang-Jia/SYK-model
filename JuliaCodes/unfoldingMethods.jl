#using Pkg
#Pkg.add("Cubature")
using Cubature
using StatsBase

function unfoldRect(rho::Function, data::Array) 
    #unfolding, using rectangles to approximate the integration of spectral density
    #rho: the spectral density functon; data: the sorted eigenvalues of a realization
    unfoldedData =  zeros(length(data))
    #the first unfolded value approximated by Integrate(rho(x), {x, -Infinity, data[1]})
    step = abs(data[2]-data[1])/10 #discretization step in the Riemann sum
    discretization =  -10*abs(data[1]) : step : data[1] #approximate -Infinity by -10*abs(data[1])
    unfoldedData[1] = sum( step*(rho.(discretization))[1:length(discretization)-1] ) 
    for i = 1: length(data)-1
        unfoldedData[i+1] = unfoldedData[i] + (data[i+1] - data[i]) * rho(data[i])
    end  
    return unfoldedData
end


function unfoldTrapz(rho::Function, data::Array) 
    #unfolding, using trapezoids to approximate the integration of spectral density
    #rho: the spectral density functon; data: the sorted eigenvalues of a realization
    unfoldedData =  zeros(length(data))
    #the first unfolded value approximated by Integrate(rho(x), {x, -Infinity, data[1]})
    step = abs(data[2]-data[1])/10 #discretization step in the Riemann sum
    discretization =  -10*abs(data[1]) : step : data[1] #approximate -Infinity by -10*abs(data[1])
    unfoldedData[1] = sum(rho.(discretization)*step) 
    for i = 1: length(data)-1
        unfoldedData[i+1] = unfoldedData[i] + (data[i+1] - data[i]) * (rho(data[i+1]) + rho(data[i]))/2
    end  
    return unfoldedData
end


function unfoldhquad(rho::Function, data::Array) #slow and not significantly better than trapezoid for SYK
    #unfolding, using trapezoids to approximate the integration of spectral density
    #rho: the spectral density functon; data: the sorted eigenvalues of a realization
    unfoldedData =  zeros(length(data))
    #the first unfolded value approximated by Integrate(rho(x), {x, -Infinity, data[1]})
     #approximate -Infinity by -10*abs(data[1])
    unfoldedData[1] = hquadrature(rho, -10*abs(data[1]), data[1],reltol=1e-7, abstol=0, maxevals=0)[1]
    for i = 1: length(data)-1
        unfoldedData[i+1] = unfoldedData[i] + hquadrature(rho, data[i], data[i+1],reltol=1e-7, abstol=0, maxevals=0)[1]
    end  
    return unfoldedData
end

function stairCase(x::Real, unfdata::Array)
    #unfdata: unfolded eigenvalues, sorted from min to max
    temp = 0 #number of eigenvalues which are smaller than x
    itr = 1
    unfdataAppended = append!(unfdata, 10^12) #makes sure itr+1 dosen't go over length(unfdata)
    while unfdataAppended[itr] < x 
        temp = temp + 1
        itr = itr +1
    end
    return temp 
end
 
function histogramData(data::Array, min::Real, max::Real, nbins::Int64) #returns  (bin's middle points, weights in bins)
    bins = range(min, stop = max, length = nbins +1)
    binMidPoints = range( (bins[1] + bins[2])/2, stop = (bins[end] + bins[end-1])/2, length = nbins  )
    histo = fit(Histogram, data, bins) 
    return (binMidPoints, histo.weights) 
end