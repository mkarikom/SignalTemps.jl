# produce the alpha matrix for cell i from equation 2
# given:
# 1) a collection of ortholog ids L
# 2) a collection of ortholog ids R
# 3) single cell data in the form of a 3d matrix
# 4)
function getalpha(df,Lig,Rec,ϵ)
    Li = filteravg(df,Lig)
    Ri = filteravg(df,Rec)
    LR = Li .* Ri'
    a = similar(LR)
    for aa in eachindex(a)
        if LR[aa] <= 1/10^ϵ
            a[aa] = 0
        else
            a[aa] = exp(-1/LR[aa])
        end
    end
    round.(a,digits=ϵ)
end

function getbeta(df,uptarg,ϵ)
    ys = filtersum(df,uptarg)
    m = length(uptarg)
    for i in eachindex(ys)
        if ys[i] <= 1/10^ϵ
            ys[i] = 0
        else
            ys[i] = exp(-m/ys[i])
        end
    end
    round.(ys,digits=ϵ)
end

function getkappa(alpha,beta,ϵ)
    k = alpha .+ beta
    for i in eachindex(k)
        if k[i] <= 1/10^ϵ
            k[i] = 0
        else
            k[i] = alpha[i] / k[i]
        end
    end
    round.(k,digits=ϵ)
end

# all up targets
# alpha * beta * kappa
function getup(df,l,r,ut,ϵ)
    # make sure the genes are expressed
    alpha = getalpha(df,l,r,ϵ)
    beta = getbeta(df,ut,ϵ)
    kappa = getkappa(alpha,beta,ϵ)
    p = alpha .* kappa .* beta
    round.(p,digits=ϵ)
end

function getupavg(df,lrt,ϵ;verbose=false)
    P = zeros(length(unique(df.barcode)),length(unique(df.barcode)))
    n = size(lrt)[1]
    for t in 1:n
        l = filterexpressed(df,lrt[t,:].ligand)
        r = filterexpressed(df,lrt[t,:].receptor)
        ut = filterexpressed(df,lrt[t,:].target)
        if all((!isnothing).([l,r,ut]))
            verbose ? println("adding triple $t of $n") : nothing
            p = getup(df,l,r,ut,ϵ)
            P += round.(p,digits=ϵ)
        else
            verbose ? println("triple $t not detected") : nothing
        end
    end

    if any(P .> 0)
        return P ./ n
    else
        verbose ? println("no lrt expression") : nothing
        return P ./ n
    end
end

function filterexpressed(df,gnames)
    nms = names(df)[findall((!isnothing).([findfirst(feat->feat==g,gnames) for g in names(df)]))]
    if length(nms) > 0
        return nms
    else
        return nothing
    end
end

# get the average
function groupavg(p,groups)
    n = size(p,1)
    n_g = length(groups)
    avg1 = Array{Float64,2}(undef,n,n_g)
    for i in 1:length(groups)
        grp = groups[i]
        avg1[:,i] = mean(p[:,grp],dims=2)
    end
    avg2 = Array{Float64,2}(undef,n_g,n_g)
    for i in 1:length(groups)
        grp = groups[i]
        avg2[i,:] = mean(avg1[grp,:],dims=1)
    end
    avg2
end

# get the groups
# find the indices of each unique value in the provided columns
# return the indices and the values
function getgroups(df,cnames...)
    nms = DataFrame([String,Vector{Int64}],[:name,:inds])
    for cn in cnames
        dfcol = df[!,cn]
        unames = unique(skipmissing(dfcol))
        for nm in unames
            nmind = []
            for i in 1:length(dfcol)
                if !ismissing(dfcol[i])
                    if dfcol[i] == nm
                        push!(nmind,i)
                    end
                end
            end
            push!(nms,(name=nm,inds=nmind))
        end
    end
    nms
end
