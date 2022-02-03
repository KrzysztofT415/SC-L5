#!/usr/bin/julia
#=
 - Author: Krzysztof Tałałaj
=#
module FileManager
export readMatrix, readVector, writeVector

using SparseArrays
using DelimitedFiles
using ...MyMap

function readMatrix(file::IOStream)
    sn, sl = split(readline(file))
    n = parse(Int64, sn)
    l = parse(Int64, sl)
    M = Dict{MyKey, Float64}()
    for line in eachline(file)
        row, column, value = split(line)
        setValue(M, parse(Int64, row), parse(Int64, column), parse(Float64, value))
    end
    M, n, l
end
export readMatrix

function readVector(file::IOStream)
    words = split(readline(file))
    n = parse(Int, words[1])
    b = Array{Float64}(undef, n)
    for (i, line) in enumerate(eachline(file))
        b[i] = parse(Float64, line)
    end
    b    
end
export readVector

function writeVector(filename::String, vector::Vector{Float64})
    writedlm(filename, vector, "\n")
end
export writeVector

end