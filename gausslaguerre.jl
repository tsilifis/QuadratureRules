"""
Gauss-Laguerre 1-dimensional rule
Author: Panos Tsilifis (pantsili@gmail.com)
Date: 09/11/2022
"""

using LinearAlgebra

function GaussLaguerre(n::Int64, odd::Bool=false)
    if odd
        n = 2^n - 1
    end
    d = [i for i=1:n-1]
    α = [2i+1 for i=0:n-1]
    Tri = LinearAlgebra.SymTridiagonal(α, d)
    ιδιοτιμες, ιδιοδιανυσματα = LinearAlgebra.eigen(Tri)
    w = ιδιοδιανυσματα[1,:] .^2
    rule = (x=ιδιοτιμες, w=w)
    return rule
end
