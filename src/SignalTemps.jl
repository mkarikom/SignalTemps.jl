module SignalTemps

using CUDA, Distributed
using LinearAlgebra
using DataFrames,CSV
using StatsBase, Distributions, Random

export getTemps

include("temps.jl") # calculate interactions

include("data.jl") # load in the data
end
