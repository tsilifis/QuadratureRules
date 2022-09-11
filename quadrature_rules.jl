"""
QuadratureRule module
Author: Panos Tsilifis (pantsili@gmail.com)
Date: 09/11/2022
"""

module QuadratureRule

    include("utils.jl")
    include("clenshawcurtis.jl")
    include("trapezoidal.jl")
    include("gausshermite.jl")
    include("gausslegendre.jl")
    include("gausslaguerre.jl")

    function QuadPoints(dim::Int64, level::Int64, sparse::Bool=true, quadrule::String="GHermite")

        @assert quadrule in ["GHermite", "GLegendre", "GLaguerre", "CC", "NC"] "Invalid rule name. Must be among 'GHermite', 'GLegendre', 'GLaguerre', 'NC' of 'CC'."
        @assert level > 0 "Level must be positive."
        @assert dim > 0 "Dimension must be positive."

        if dim==1
            if quadrule == "GHermite"
                grid = GaussHermite(level, sparse)
                x, w = grid.x, grid.w
            elseif quadrule == "GLegendre"
                grid = GaussLegendre(level, sparse)
                x, w = grid.x, grid.w
            elseif quadrule == "GLaguerre"
                grid = GaussLaguerre(level, sparse)
                x, w = grid.x, grid.w
            elseif quadrule == "NC"
                grid = Trapezoidal(level)
                x, w = grid.x, grid.w
            else
                grid = ClenshawCurtis(level)
                x, w = grid.x, grid.w
            end
            return x, w
        else
            d0 = dim - 1
            if quadrule == "GHermite"
                Q1 = [GaussHermite(i, sparse) for i=1:level+d0-1]
                H_rule = GaussHermite(1, sparse)
                D1 = [[H_rule.x, H_rule.w]]
            elseif quadrule == "GLegendre"
                Q1 = [GaussLegendre(i, sparse) for i=1:level+d0-1]
                L_rule = GaussLegendre(1, sparse)
                D1 = [[L_rule.x, L_rule.w]]
            elseif quadrule == "GLaguerre"
                Q1 = [GaussLaguerre(i, sparse) for i=1:level+d0-1]
                L_rule = GaussLaguerre(1, sparse)
                D1 = [[L_rule.x, L_rule.w]]
            elseif quadrule == "CC"
                Q1 = [ClenshawCurtis(i) for i=1:level+d0-1]
                CC_rule = ClenshawCurtis(1)
                D1 = [[CC_rule.x, CC_rule.w]]
            else
                Q1 = [Trapezoidal(i) for i=1:level+d0-1]
                NC_rule = Trapezoidal(1)
                D1 = [[NC_rule.x, NC_rule.w]]
            end
            for i=2:length(Q1)
                D1 = hcat(D1, RuleDiff(Q1[i], Q1[i-1]))
            end

            MI = MultiIndex(level+d0-1, d0)
            THETA = zeros(dim)
            W = zeros(1)
            for i=1:length(MI[:, 1])
                mi = MI[i, :]
                k_abs = sum(mi)
                x = [y for y in Q1[level+d0-k_abs].x]
                w = [y for y in Q1[level+d0-k_abs].w]
                for j=1:d0
                    x = [vcat(u, v) for u in x for v in D1[mi[j]][1]]
                    w = [hcat(u, v) for u in w for v in D1[mi[j]][2]]
                end
                for i=1:length(x)
                    THETA = hcat(THETA, x[i])
                    W = hcat(W, prod(w[i]))
                end
            end
            if level>1
                W = reshape(W[1, 2:end], 1, length(W)-1)
                THETA = THETA[:, 2:end]'
                theta_list = sort([THETA[i, :] for i=1:length(THETA[:, 1])])
                #theta_uni = unique(THETA, dims=1)
                theta_uni_list = unique(theta_list)
                theta_uni = mapreduce(permutedims, vcat, theta_uni_list)
                ind = unique(i->theta_list[i], 1:length(theta_list))
                c = [count(==(i), theta_list) for i in unique(theta_list)]
                inv = [findall(theta_uni_list->theta_uni_list==i, theta_uni_list) for i in theta_list]

                #println("theta uni")
                #println(theta_uni)
                #println("ind")
                #println(ind)
                #println("inv")
                #println(inv)
                #println("c")
                #println(c)
                locs = [j for j=1:length(c) if c[j] > 1]
                w_uni = [W[1, i] for i in ind]
                for j in locs
                    loc_j = [k for k=1:length(inv) if inv[k][1]==j]
                    if length(loc_j) > 0
                        w_uni[j] = sum([W[1, locj] for locj in loc_j])
                    #else
                    #    w_uni[j] = 0.
                    end
                end
                locs_0 = findall(theta_uni->abs(theta_uni)<1e-16, theta_uni)
                for loc in locs_0
                    theta_uni[loc] = 0.
                end
                x_final = theta_uni
                w_final = w_uni
            else
                x_final = THETA[:, 2:end]'
                w_final = reshape(W[1, 2:end], 1, length(W)-1)
            end
            return x_final, w_final
        end
    end
end
