# return the subset of expressed markers
function filterdf(df,colnames)
    filtnames = []
    for cn in colnames
        if any(occursin.(cn,names(df)))
            push!(filtnames,cn)
        end
    end
    filtnames
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
