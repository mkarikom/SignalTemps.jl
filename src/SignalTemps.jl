module SignalTemps

using CUDA, Distributed
using LinearAlgebra
using DataFrames,Query
using StatsBase, Distributions, Random
using UMAP
using LightGraphs,MetaGraphs,PCquery

export getalpha, getbeta, getkappa, gettemps, safediv, getup, getupavg, groupavg, getgroups
export filterdf,filteravg,filtersum,umapexpr,sublabels,fplot,filterexpressed
export enumerateLRT,lrtDist,lrtManual

include("temps.jl") # calculate interactions

include("data.jl") # load in the data

include("lrtGraph.jl") # pcquery graph based lrt ident

include("lrtManual.jl") # manual lrt ident, uses OrthoDB ortholog (eg danio) search from PCquery given primary (eg human) pathway genes
end
