#!/usr/bin/julia
#=
 - Author: Krzysztof Tałałaj
=#
module blocksys

using ...MyMap

function calculateRightSide(M::Dict{MyKey, Float64}, n::Int64, l::Int64)
    b = zeros(Float64, n)
    for i in 1:n, j in max(1, i - (2 + l)):min(n, i + l)
        b[i] += getValue(M, i, j)
    end
    b
end
export calculateRightSide


function gauss!(M::Dict{MyKey, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    for k in 1:n-1
        for i in k+1:min(n, k + l + 1)
            z = getValue(M, i, k) / getValue(M, k, k)
            setValue(M, i, k, 0.0)

            for col in k+1:min(n, k + l + 1)
                setValue(M, i, col, getValue(M, i, col) - z * getValue(M, k, col))
            end

            b[i] -= z * b[k]
        end
    end
end

function gaussPivoted!(M::Dict{MyKey, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    mains = collect(1:n)

    for k in 1:n-1
        max_v, max_i = 0, 0

        for i in k:min(n, k + l + 1)
            if abs(getValue(M, mains[i], k)) > max_v
                max_v = abs(getValue(M, mains[i], k))
                max_i = i
            end
        end

        mains[max_i], mains[k] = mains[k], mains[max_i]

        for i in k+1:min(n, k + l + 1)
            z = getValue(M, mains[i], k) / getValue(M, mains[k], k)
            setValue(M, mains[i], k, 0.0)

            for j in k+1:min(n, k + 2 * l)
                setValue(M, mains[i], j, getValue(M, mains[i], j) - z * getValue(M, mains[k], j))
            end
            b[mains[i]] -= z * b[mains[k]]
        end
    end

    return mains
end

function gaussLU!(M::Dict{MyKey, Float64}, n::Int64, l::Int64)
    L = Dict{MyKey, Float64}()
    for k in 1:n-1
        setValue(L, k, k, 1.0)
        for i in k+1:min(n, k + l + 1)
            z = getValue(M, i, k) / getValue(M, k, k)
            setValue(L, i, k, z)
            setValue(M, i, k, 0.0)
            for j in k+1:min(n, k + l)
                setValue(M, i, j, getValue(M, i, j) - z * getValue(M, k, j))
            end
        end
    end
    setValue(L, n, n, 1.0)
    L
end

function gaussPivotedLU!(M::Dict{MyKey, Float64}, n::Int64, l::Int64)
    L, mains = Dict{MyKey, Float64}(), collect(1:n)

    for k in 1:n-1
        max_v, max_i = 0, 0

        for i in k:min(n, k + l + 1)
            if abs(getValue(M, mains[i], k)) > max_v
                max_v = abs(getValue(M, mains[i], k))
                max_i = i
            end
        end

        mains[max_i], mains[k] = mains[k], mains[max_i]

        for i in k+1:min(n, k + l + 1)
            z = getValue(M, mains[i], k) / getValue(M, mains[k], k)

            setValue(L, mains[i], k, z)
            setValue(M, mains[i], k, 0.0)

            for j in k+1:min(n, k + 2 * l)
                setValue(M, mains[i], j, getValue(M, mains[i], j) - z * getValue(M, mains[k], j))
            end
        end
    end

    return L, mains
end

function runGauss(M::Dict{MyKey, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    gauss!(M, b, n, l)
    x_n = zeros(Float64, n)

    for i in n:-1:1
        x_i = 0
        for j in i+1:min(n, i + l)
            x_i += getValue(M, i, j) * x_n[j]
        end
        x_n[i] = (b[i] - x_i) / getValue(M, i, i)
    end
    return x_n
end
export runGauss

function runPivotedGauss(M::Dict{MyKey, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    mains = gaussPivoted!(M, b, n, l)
    x_n = zeros(Float64, n)

    for k in 1:n-1
        for i in k+1:min(n, k + 2 * l)
            b[mains[i]] -= getValue(M, mains[i], k) * b[mains[k]]
        end
    end

    for i in n:-1:1
        x_i = 0
        for j in i+1:min(n, i + 2 * l)
            x_i += getValue(M, mains[i], j) * x_n[j]
        end
        x_n[i] = (b[mains[i]] - x_i) / getValue(M, mains[i], i)
    end

    return x_n
end
export runPivotedGauss

function runGaussLU(M::Dict{MyKey, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    L = gaussLU!(M, n, l)
    x_n = zeros(Float64, n)

    for k in 1:n-1
        for i in k+1:min(n, k + l + 1)
            b[i] -= getValue(L, i, k) * b[k]
        end
    end

    for i in n:-1:1
        x_i = 0
        for j in i+1:min(n, i + l)
            x_i += getValue(M, i, j) * x_n[j]
        end
        x_n[i] = (b[i] - x_i) / getValue(M, i, i)
    end

    return x_n
end
export runGaussLU

function runPivotedGaussLU(M::Dict{MyKey, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    L, mains = gaussPivotedLU!(M, n, l)
    x_n = zeros(Float64, n)

    for k in 1:n-1
        for i in k+1:min(n, k + l + 1)
            b[mains[i]] -= getValue(L, mains[i], k) * b[mains[k]]
        end
    end

    for i in n:-1:1
        x_i = 0
        for j in i+1:min(n, i + 2 * l)
            x_i += getValue(M, mains[i], j) * x_n[j]
        end
        x_n[i] = (b[mains[i]] - x_i) / getValue(M, mains[i], i)
    end

    return x_n
end
export runPivotedGaussLU

end
