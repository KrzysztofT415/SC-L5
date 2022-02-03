#!/usr/bin/julia
#=
 - Author: Krzysztof Tałałaj
=#
println("Importing packages")

include("src/myMap.jl")
using .MyMap
include("src/matrixgen.jl")
include("src/blocksys.jl")
using .blocksys
using .matrixgen

using Test
# using Plots
using LinearAlgebra

function testSolver(solver, solverName::String, l::Int64, a::Int64, b::Int64, d::Int64, i::Int64)
    memory::Vector{Float64}, time::Vector{Float64}, r::Float64 = [], [], 0
    for n in a:d:b
        t, m = 0, 0
        for _ in 1:i
            M = blockmat(n, l, 1.0)
            b = calculateRightSide(M, n, l)
            x, dt, dm = @timed solver(M, b, n, l)
            @test isapprox(x, ones(Float64, n))
            r += norm(ones(Float64, n) - x) / norm(x)
            t += dt
            m += dm
        end
        println(solverName, "\t", l, " ", n, "\t", m / i / 1024, "  \t", t / i)
        push!(time, t / i)
        push!(memory, m / i / 1024)
    end
    time, memory, r / i / length(memory)
end

println("Starting main")
println("solver          l n\tmemory\t\ttime")

a, b, d, i = 10000, 10000, 2, 3
xAxis = [n for n in a:d:b]
R = Vector{Float64}()

t, m, r = testSolver(runGauss, "Gauss         ", 8, a, b, d, i)
# p1 = plot(xAxis, t, label = "Gauss", xlims = (a, b), lw = 2, title="time l=8")
# p2 = plot(xAxis, m, label = "Gauss", xlims = (a, b), lw = 2, title="memory l=8")
push!(R, r)

t, m, r = testSolver(runPivotedGauss, "PivotedGauss  ", 8, a, b, d, i)
# plot!(p1, xAxis, t, label = "PivotedGauss", xlims = (a, b), lw = 2)
# plot!(p2, xAxis, m, label = "PivotedGauss", xlims = (a, b), lw = 2)
push!(R, r)

t, m, r = testSolver(runGaussLU, "GaussLU       ", 8, a, b, d, i)
# plot!(p1, xAxis, t, label = "GaussLU", xlims = (a, b), lw = 2)
# plot!(p2, xAxis, m, label = "GaussLU", xlims = (a, b), lw = 2)
push!(R, r)

t, m, r = testSolver(runPivotedGaussLU, "PivotedGaussLU", 8, a, b, d, i)
# plot!(p1, xAxis, t, label = "PivotedGaussLU", xlims = (a, b), lw = 2)
# plot!(p2, xAxis, m, label = "PivotedGaussLU", xlims = (a, b), lw = 2)
push!(R, r)


t, m, r = testSolver(runGauss, "Gauss         ", 8, 100000, 100000, d, i)
# p3 = plot(xAxis, t, label = "Gauss", xlims = (a, b), lw = 2, title="time l=8")
# p4 = plot(xAxis, m, label = "Gauss", xlims = (a, b), lw = 2, title="memory l=8")
push!(R, r)

println(R)
# plot(p1, p2, p3, p4, layout = (2, 2), size = (1920, 1080))