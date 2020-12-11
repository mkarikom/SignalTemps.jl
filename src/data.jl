# return the subset of expressed markers
function filterdf(df,colnames)
    filtnames = []
    for cn in colnames
        ff = findfirst(n->n==cn,names(df))
        if !isnothing(ff)
            push!(filtnames,cn)
        end
    end
    unique(filtnames)
end

# get the average over a subset of markers
function filteravg(df,colnames)
    cn = filterdf(df,colnames)
    mean(eachcol(df[!,Symbol.(cn)]))
end

# get the sum over a subset of markers
function filtersum(df,colnames)
    cn = filterdf(df,colnames)
    sum(eachcol(df[!,Symbol.(cn)]))
end

# nms is feature columns for umap, all others interpreted as labels during stack
function umapexpr(df,nms,params)
    inds = [all(map(n->n!=nm,nms)) for nm in names(df)]
    um = umap(Matrix(df[!,inds])',
                n_neighbors=params[:nneigh],
                min_dist=params[:mindist],
                n_epochs=params[:nepochs])

    # scatter(res_jl[1,:], res_jl[2,:], zcolor=mnist_y,
    #         title="MNIST: Julia UMAP", marker=(2, 2, :auto, stroke(0)))
    #         inds = [all(map(n->n!=nm,nms)) for nm in names(df)]

end


# feature plot with dimensionality reduction
# df is a dataframe with shape samples x features
# ids are identifiers, such that names(df) \ ids  = the features we care about
# trans is a transformation of the samples x features passed to the dimred call
# dimred is a function that outputs a 2d feature vector
# the feature vector is added to df as feat_x and feat_y
function fplot(df,ids,dimred,trans,args)
    # NOT IMPLEMENTED
    dat = df[!,[all(map(n->n!=nm,ids)) for nm in names(df)]]
    reduced = eval(Expr(:call,dimred,trans(dat),args))
end

# given a list of at least 2 vectors, iteratively replace nonmissing values with right associativity
function sublabels(labels...)
    @assert length(labels) >= 2 "less than 2 label vectors provided"
    lab = labels[1]
    for l in 2:length(labels)
        for ll in 1:length(lab)
            if !ismissing(labels[l][ll])
                local lab[ll] = labels[l][ll]
            end
        end
    end
    lab
end
