function lrtManual(upids,orthos,dist)
	df = DataFrame(fill(Vector{String},3),Symbol.(keys(upids)))
	# df = []
	for l in upids[:ligand]
		for r in upids[:receptor]
			for t in upids[:target]
				ligs = filter(row->row.entId==l&&row.orthoDist==string(dist),orthos)
				recs = filter(row->row.entId==r&&row.orthoDist==string(dist),orthos)
				targs = filter(row->row.entId==t&&row.orthoDist==string(dist),orthos)
				if size(ligs,1) > 0 && size(recs,1) > 0 && size(targs,1) > 0
					newrow = Dict(:ligand=>collect(skipmissing(ligs.orthoHGNC)),
						 		 :receptor=>collect(skipmissing(recs.orthoHGNC)),
								 :target=>collect(skipmissing(targs.orthoHGNC)))
				    push!(df,newrow)
				end
			end
		end
	end
	df
end
