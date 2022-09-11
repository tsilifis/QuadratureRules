"""
Newton-Cotes 1-dimensional rule
Author: Panos Tsilifis (pantsili@gmail.com)
Date: 09/11/2022
"""

function Trapezoidal(l::Int64)
    if l==1
        x = [0.]
        w = [1.]
    else
        n = 2^(l-1) + 1
        x = [i*2^(2. -l)-1. for i=0:n-1]
        w = [2^(1. -l) for i=0:n-1]
        w[begin] /= 2.
        w[end] /= 2.
    end
    rule = (x=x, w=w)
    return rule
end
