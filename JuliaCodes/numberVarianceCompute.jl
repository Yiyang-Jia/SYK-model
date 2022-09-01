function  intervalGeneration(dataDim::Int64 ,  intervalLength::Real , overlapFrac::Real) 
    # dataDim : length of the  input column vector
    # the function creats overlaping intervals spanning [0, dataDim], i.e., this assumes the spectrum is properly unfolded so that the average position of the first eigenvalue is near 0
    if overlapFrac > 0.5
      return "Overlap fraction must be smaller than 0.5."
    else
      overlapLen = intervalLength * overlapFrac
      nIntervals =  floor(Int, (dataDim-overlapLen)/( ((1 - overlapFrac)*intervalLength) ) ) 
      intervals = Array{Tuple{Real, Real}}(UndefInitializer(), nIntervals)
      for k = 1 : nIntervals 
          intervals[k] = ( (k-1) * (1-overlapFrac) * intervalLength , (k-1) * (1-overlapFrac) * intervalLength + intervalLength )
      end
      return  intervals
    end
end

  function nVcompute(data::Array{Float64,1} , intervalLength::Real , overlapFrac::Real) #compute the number variance of data
     #overlapFrac: fraction of overlaps between neighboring intervals, should be smaller than 0.5
     #data: the sorted spectrum of one single realization (a column vector), with average spacing unfolded to 1. Smallest eigenvalues should be close to 0 (use unfoldingMethods.jl)
     #length of interval, equal to average number of levels when the data is unfolded to have unit average spacing
     if overlapFrac > 0.5
       return "Overlap fraction must be smaller than 0.5."
     else
        dataDim = length(data)
        intervals = intervalGeneration(dataDim, intervalLength, overlapFrac)
        nIntervals = length(intervals) #count the number of intervals
        numbersInIntervals = zeros(Int64, nIntervals)
        for k = 1 : dataDim
          estimator = data[k]/ ( (1-overlapFrac) * intervalLength )
          estimatePos = floor(Int, estimator)
          if  0 < estimatePos < nIntervals &&  estimator < estimatePos + overlapFrac/( 1- overlapFrac )  #if the data point is in an overlapped region
            numbersInIntervals[estimatePos] =  numbersInIntervals[estimatePos] +1
            numbersInIntervals[estimatePos+1] =  numbersInIntervals[estimatePos+1] +1
          elseif -1 < estimatePos < nIntervals  &&  estimator > estimatePos + overlapFrac/( 1- overlapFrac )                                                           #if the data point is in an un-overlapped region
            numbersInIntervals[estimatePos+1] =  numbersInIntervals[estimatePos+1] +1
          end
        end
     end
    #return numbersInIntervals   # test translational invariance
    return sum(numbersInIntervals)/nIntervals
    #return sum(numbersInIntervals.^2)/nIntervals - intervalLength^2
  end

  function nVcomputeNoOverlap(data::Array{Float64,1} , intervalLength::Real) #compute the number variance of data
    #overlapFrac: fraction of overlaps between neighboring intervals, should be smaller than 0.5
    #data: the sorted spectrum of one single realization (a column vector), with average spacing unfolded to 1. Smallest eigenvalues should be close to 0 (use unfoldingMethods.jl)
    #length of interval, equal to average number of levels when the data is unfolded to have unit average spacing
  
    overlapFrac = 0.0
    dataDim = length(data)
    intervals = intervalGeneration(dataDim, intervalLength, overlapFrac)
    nIntervals = length(intervals) #count the number of intervals
    numbersInIntervals = zeros(Int64, nIntervals)
    for k = 1 : dataDim
      estimator = data[k]/ ( (1-overlapFrac)*intervalLength )
      estimatePos = floor(Int, estimator)
      if -1 < estimatePos < nIntervals
      numbersInIntervals[estimatePos+1] =  numbersInIntervals[estimatePos+1] + 1
      end
    end
  
   #return numbersInIntervals   # test translational invariance
   return sum( (numbersInIntervals - intervalLength * ones(nIntervals))[3:end-2].^2 )/(nIntervals-4) #exclude the first few and the last few intervals
   #return (sum(numbersInIntervals[2:end-1])/(nIntervals-2) - intervalLength)^2  #exclude the first and the last interval
 end