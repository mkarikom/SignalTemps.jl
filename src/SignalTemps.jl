module SignalTemps

using CUDA, Distributed
using LinearAlgebra
using DataFrames,Query
using StatsBase, Distributions, Random

export getalpha, gettemps, gettargs
export filterdf,filteravg

include("temps.jl") # calculate interactions

include("data.jl") # load in the data
end
