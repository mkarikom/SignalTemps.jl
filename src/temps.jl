# produce the alpha matrix for cell i from equation 2
# given:
# 1) a collection of ortholog ids L
# 2) a collection of ortholog ids R
# 3) single cell data in the form of a 3d matrix
# 4)
function getalpha(df,Lig,Rec)
    Li = filteravg(df,Lig)
    Ri = filteravg(df,Rec)
    alpha = exp.(-1 ./ (Li .* Ri'))
end

function getbeta(df,uptarg)
    ys = filtersum(df,uptarg)
    m = size(ys)[2]
    beta = exp.(-m ./ ys)
end


# repeat Li over rows
function gettemps(L,n)

end

function gettargs(row)
    for r in row
        println("r")
    end
end
