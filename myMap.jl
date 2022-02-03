#!/usr/bin/julia
#=
 - Author: Krzysztof Tałałaj
=#
module MyMap

struct MyKey
    x::Int
    y::Int
end
export MyKey

function Base.isequal(A::MyKey, B::MyKey)
    A.x == B.x && A.y == B.y
end

function Base.hash(A::MyKey)
    hash(A.x + A.y)
end

function getValue(A::Dict{MyKey, Float64}, x::Int64, y::Int64)
    if (haskey(A, MyKey(x, y)))
        return A[MyKey(x, y)]
    end
    return 0
end
export getValue

function setValue(A::Dict{MyKey, Float64}, x::Int64, y::Int64, v::Float64)
    A[MyKey(x, y)] = v
end
export setValue

end