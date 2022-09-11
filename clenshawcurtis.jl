"""
Clenshaw-Curtis 1-dimensional rule
Author: Panos Tsilifis (pantsili@gmail.com)
Date: 09/11/2022
"""

function ClenshawCurtis(l::Int64)
    if l == 1
        x = [-cos(π/2)]
        w = [2.]
    else
        n = 2^(l-1) + 1
        x = [-cos(π*i/(n-1.)) for i=0:n-1]
        w = zeros(n)
        w[begin] = w[end] = 1. / (n*(n-2.))
        for i=2:n-1
            s = [cos(2π*(i-1)*(j+1)/(n-1.)) / (1. - 4*(j+1.)^2) for j=0:Int(floor((n-1)/2))-1]
            s[end] /= 2.
            w[i] = 2 * (1. + 2*sum(s)) / (n-1.)
        end
    end
    rule = (x=x, w=w/2)
    return rule
end
