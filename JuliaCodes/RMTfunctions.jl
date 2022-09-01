using SpecialFunctions
## Spectral form factors, with average level spacing being 1 
function sffGOE(t::Real) #t>0
    x = t / (2* pi)
    if x < 0
        return "Time must be postive"
    elseif x < 1
        return 2*x - x*log(1 + 2*x)
    else
        return 2 - x * log( (2*x + 1)/(2*x - 1) )
    end
end

function sffGUE(t::Real) #t>0
    x = t / (2* pi)
    if x < 0
        return "Time must be postive"
    elseif x < 1
        return x
    else
        return 1
    end
end

function sffGSE(t::Real) #t>0
    x = t / (2* pi)
    if x < 0
        return "Time must be postive"
    elseif x < 2
        return x/2 - x/4 * log( abs(1 - x) )
    else
        return 1
    end
end
sffplots=plot(0:0.01:11 , [sffGOE, sffGUE, sffGSE], title = "RMT spectral form factors", label=["GOE" "GUE" "GSE"], linewidth=2, legend=:bottomright)
savefig(sffplots, "/Users/yiyang/Desktop/SYKcodes/JuliaCodes/sffplots.png")

######################################################################################################################################################################################################################################################################################################
## number variance

function nVGOE(s::Real)
    if s == 0
        return 0
    else
        eulerGamma = 0.5772156649015329
        firstTerm = 2/pi^2 * ( log( 2*pi*s ) + eulerGamma - pi^2/8 - cosint( 2*pi*s ) )
        secondTerm =  4*s/pi * ( ( 1 + pi^2*s -cos(2*pi*s) )/(2*pi*s) - sinint(2*pi*s) )
        thridTerm = 1/pi^2 * ( ( pi - 2 * sinint(pi*s) )/2 )^2
        return  firstTerm + secondTerm + thridTerm
    end
end

function nVGUE(s::Real)
    if s == 0
        return 0
    else
        eulerGamma = 0.5772156649015329
        firstTerm = 1/pi^2 * ( log( 2*pi*s ) + eulerGamma - cosint( 2*pi*s ) )
        secondTerm =  2*s/pi * ( ( 1 + pi^2*s -cos(2*pi*s) )/(2*pi*s) - sinint(2*pi*s) )
        return  firstTerm + secondTerm
    end
end

function nVGSE(s::Real)
    if s == 0
        return 0
    else
        eulerGamma = 0.5772156649015329
        firstTerm = 1/(2*pi^2) * ( log( 4*pi*s ) + eulerGamma - cosint( 4*pi*s ) )
        secondTerm =  2*s/pi * ( ( 1 + pi^2*2*s -cos(4*pi*s) )/(4*pi*s) - sinint(4*pi*s) )
        thridTerm = 1/(4*pi^2) * ( sinint(2*pi*s ) )^2
        return  firstTerm + secondTerm + thridTerm
    end
end
nvplots= plot(0:0.01:6 , [nVGOE, nVGUE, nVGSE], title = "RMT number variances", label=["GOE" "GUE" "GSE"],linewidth=2, legend=:bottomright)
savefig(nvplots, "/Users/yiyang/Desktop/SYKcodes/JuliaCodes/nvplots.png")

######################################################################################################################################################################################################################################################################################################
## cluster functions (= - (DensityDensityConnectedCorrelation - deltaFunction)) in the microscopic limit (namely after unfolding, and don't look too long-range)
function twoPointClusterGOE(r::Real)
    if r == 0
        return 1
    else
        return   ( 1/2 - sinint( pi * r )/pi ) * ( cos( pi* r)/r - sin( pi* r )/( pi* r^2) ) + ( sin(pi*r)/(pi*r) )^2

    end
end

function twoPointClusterGUE(r::Real)
    if r == 0
        return 1
    else
        return  ( sin(pi*r)/(pi*r) )^2
    end
end

function twoPointClusterGSE(r::Real)
    if r == 0
        return 1
    else
        return  ( sin(2*pi*r)/(2*pi*r) )^2  - ( cos(2*pi*r)/(2*pi*r) - sin(2*pi*r)/(4*pi^2*r^2) ) * sinint(2*pi*r)
    end
end
#plot(0:0.01:6 , [twoPointClusterGOE, twoPointClusterGUE, twoPointClusterGSE], title = "Two-point cluster functions", label=["GOE" "GUE" "GSE"])
######################################################################################################################################################################################################################################################################################################