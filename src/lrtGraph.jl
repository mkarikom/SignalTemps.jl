## find all possible triples based on target paths from sParam[:d0]
function enumerateLRT(g,paths,drange;verbose=false)
    up0 = []
    ens0 = []
    hgnc0 = []
    upOrth = []
    hgncOrth = []
    pathInd = []
    count = 0
    for p in 1:length(paths)
        count+=1
        verbose ? println("running path $count") : nothing
        # sleep(0.5)
        v0 = paths[p][:v0]
        for v1 in paths[p][:v1]
            pv0 = props(g,v0)
            pv1 = props(g,v1)
            ptarg = props(g,paths[p][:targProt])
            push!(up0,Dict(
                    Symbol(pv0[:roleLR])=>pv0[:entId],
                    Symbol(pv1[:roleLR])=>pv1[:entId],
                    :target=>ptarg[:entId]))
            push!(ens0,Dict(
                    Symbol(pv0[:roleLR])=>pv0[:ensId],
                    Symbol(pv1[:roleLR])=>pv1[:ensId],
                    :target=>ptarg[:ensId]))
            push!(hgnc0,Dict(
                    Symbol(pv0[:roleLR])=>pv0[:gname],
                    Symbol(pv1[:roleLR])=>pv1[:gname],
                    :target=>ptarg[:gname]))

            orthv0hgnc = Dict()
            orthv0upid = Dict()
            for d in drange
                # conditional to avoid error due to bad jld2 save
                # if typeof(props(g,v0)[:orthoDist][d][:members]) <: NTuple && typeof(props(g,v1)[:orthoDist][d][:members]) <: NTuple && typeof(props(g,p[:targProt])[:orthoDist][d][:members]) <: NTuple
                    odv0 = getOrthDist(g,v0,d)
                    odv1 = getOrthDist(g,v1,d)
                    odt = getOrthDist(g,paths[p][:targProt],d)
                    push!(orthv0hgnc,d=>Dict(Symbol(pv0[:roleLR])=>odv0.hgnc,
                                             Symbol(pv1[:roleLR])=>odv1.hgnc,
                                             :target=>odt.hgnc))
                    push!(orthv0upid,d=>Dict(Symbol(pv0[:roleLR])=>odv0.upid,
                                             Symbol(pv1[:roleLR])=>odv1.upid,
                                             :target=>odt.upid))
                # end
            end
            push!(hgncOrth,orthv0hgnc)
            push!(upOrth,orthv0upid)
            push!(pathInd,p)
        end
    end
    (up0=unique(up0),ens0=unique(ens0),hgnc0=unique(hgnc0),hgncOrth=hgncOrth,upOrth=upOrth,pathInd=pathInd)
end

# create a triples dataframe from output of enumerateLRT
function lrtDist(arr,dist,marker)
    df = DataFrame(ligand=Vector{Vector{String}}(),receptor=Vector{Vector{String}}(),target=Vector{Vector{String}}(),pathInd=Vector{Int64}())
    for a in 1:length(arr[marker])
        push!(df,Dict(arr[marker][a][dist]...,:pathInd=>arr[:pathInd][a]))
    end
    unique(df[!,[:ligand,:receptor,:target]])
end
