function groundStateE(sig, eta) #emax
    return sqrt(4*sig^2/(1 - eta))
end 

function suppressionFactor(n, q) #eta for syk
    temp = 0
    for r = 0:q
        temp = temp + (-1)^(q + r) * binomial(q, r) * binomial(n - q, q - r)
    end
    return temp/binomial(n,q)
end 

function qHdensityNormalized(x::Real, eta::Real, sig::Number, ntrunc::Int64) 
    emax = groundStateE(sig, eta)
    productFactor = 1 
    normFactor = 1/(pi*sig) * (1+eta) * sqrt( 1 - eta )
    for n = 1: ntrunc
        productFactor = productFactor * ( 1- 4*x^2/( emax^2 * (eta^n + eta^(-n)+2) ) )
        normFactor = normFactor *  ( 1 - eta^(2*n + 2) ) / ( 1 - eta^(2*n + 1) )
    end
    if abs(x) < emax
        return normFactor * sqrt(1 - x^2/emax^2)* productFactor
    else
        return 0
    end
end

function  qHpolynomial(x::Real,  n::Int64, eta::Real, sig::Number)
    y = x/sig
    function h(y, eta, n)
        if n == 0
            return 1
        elseif n == 1
            return y 
        else
            firstTerm = y * h(y, eta, n-1) 
            prefactorOfSecondTerm = 0
            for i = 0 : n-2
                prefactorOfSecondTerm = prefactorOfSecondTerm + eta^i
            end
            secondTerm = - prefactorOfSecondTerm * h(y, eta, n-2)
            return firstTerm + secondTerm
        end
    end
    return h(y, eta, n)

    
end
# test validity:
# x= 0.21
# eta = 0.88
# (x^2 - 1, qHpolynomial(x, 2, eta, 1))

# (-(1 + eta)* x + x *(-1 + x^2),  qHpolynomial(x,3, eta, 1))
