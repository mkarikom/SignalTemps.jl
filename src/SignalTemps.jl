module SignalTemps

using CUDA, Distributed
using LinearAlgebra
using DataFrames,Query
using StatsBase, Distributions, Random
using UMAP
using LightGraphs,MetaGraphs

export getalpha, getbeta, getkappa, gettemps, safediv, getup, getupavg, groupavg, getgroups
export filterdf,filteravg,filtersum,umapexpr,sublabels,fplot,filterexpressed
export lrtManual

include("temps.jl") # calculate interactions

include("data.jl") # load in the data

include("lrtManual.jl") # manual lrt ident, uses OrthoDB ortholog (eg danio) search from PCquery given primary (eg human) pathway genes
end
