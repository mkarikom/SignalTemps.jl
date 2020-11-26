module SignalTemps

using CUDA, Distributed
using LinearAlgebra
using DataFrames,Query
using StatsBase, Distributions, Random

export getalpha, getbeta, getkappa, gettemps, safediv, getup, getupavg
export filterdf,filteravg,filtersum

include("temps.jl") # calculate interactions

include("data.jl") # load in the data
end
