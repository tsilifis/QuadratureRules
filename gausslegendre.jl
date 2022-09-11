"""
Gauss-Legendre 1-dimensional rule
Author: Panos Tsilifis (pantsili@gmail.com)
Date: 09/11/2022
"""

using LinearAlgebra

function GaussLegendre(n::Int64, odd::Bool=false)
    if odd
        n = 2^n - 1
    end
    d = sqrt.([i^2/((2i+1.)*(2i-1.)) for i=1:n-1])
    Tri = LinearAlgebra.SymTridiagonal(zeros(n), d)
    ιδιοτιμες, ιδιοδιανυσματα = LinearAlgebra.eigen(Tri)
    μ = length(ιδιοτιμες)
    w = ιδιοδιανυσματα[1,:] .^2
    if (μ-1)%2 == 0
        ιδιοτιμες[Int(ceil(μ/2))] = 0.
    end
    rule = (x=ιδιοτιμες, w=w)
    return rule
end
