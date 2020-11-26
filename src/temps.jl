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
    alpha = getalpha(df,l,r,ϵ)
    beta = getbeta(df,ut,ϵ)
    kappa = getkappa(alpha,beta,ϵ)
    p = alpha .* kappa .* beta
    round.(p,digits=ϵ)
end

function getupavg(df,lrt,ϵ)
    P = nothing
    n = size(lrt)[1]
    for t in 1:n
        println("running triple $t of $n")
        l = lrt[t,:].ligand
        r = lrt[t,:].receptor
        ut = lrt[t,:].target
        if isnothing(P)
            P = getup(df,l,r,ut,ϵ)
        else
            P += getup(df,l,r,ut,ϵ)
        end
    end
    round.(P ./ n, digits=ϵ)
end
