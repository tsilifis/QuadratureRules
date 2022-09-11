"""
Gauss-Hermite 1-dimensional rule
Author: Panos Tsilifis (pantsili@gmail.com)
Date: 09/11/2022
"""

using LinearAlgebra

function GaussHermite(n::Int64, odd::Bool=false)
    if odd
        n = 2^n - 1
    end
    d = sqrt.(range(1, n-1))
    Tri = LinearAlgebra.SymTridiagonal(zeros(n), d)
    ιδιοτιμες, ιδιοδιανυσματα = LinearAlgebra.eigen(Tri)
    if odd
        ιδιοτιμες[Int(ceil(n/2))] = 0.
    end
    w = ιδιοδιανυσματα[1,:] .^2
    rule = (x=ιδιοτιμες, w=w)
    return rule
end
