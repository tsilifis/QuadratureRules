"""
Util functions for multi-dimensional quadratures. 
Author: Panos Tsilifis (pantsili@gmail.com)
Date: 09/11/2022
"""

function nchoosek(n::Int64, k::Int64)
    return Int(factorial(big(n)) / (factorial(big(k))*factorial(big(n-k))))
end

function MultiIndex(order::Int64, dim::Int64)
    @assert dim > 0 "Expected positive dimension"
    q_num = [nchoosek(dim+i-1, i) for i=0:order]
    mul_ind = zeros(Int64, 1, dim)
    mul_ind = [mul_ind; Matrix(1I, dim, dim)]
    eye = Matrix(1I, dim, dim)
    ind = [1 for i=1:dim]
    terms = []
    for j=1:order-1
        ind_new = zeros(Int64, 0)
        for i=0:dim-1
            a0 = eye[sum(ind[begin:i])+1:end, :]
            a0[:, i+1] .+= 1
            mul_ind = vcat(mul_ind, a0)
            ind_new = vcat(ind_new, length(a0[:, 1]))
        end
        ind = ind_new
        eye = mul_ind[sum(q_num[begin:j+1])+1:end, :]
    end
    terms = [mul_ind[i, :] for i=1:length(mul_ind[:, 1]) if minimum(mul_ind[i, :])>0]
    terms = hcat(terms...)'
    return terms
end

function RuleDiff(rule1, rule2)
    return [rule1.x; rule2.x], [rule1.w; -rule2.w]
end
