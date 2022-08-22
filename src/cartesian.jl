
"""
Sample each point once, spread uniform randomly over the dynamic axis
"""
function uniform(shape::NTuple{N, Integer}) where N
	indices = vec(collect(CartesianIndices(shape)))
	return shuffle!(indices)
end

function poisson_disk()
	error("Not implemented")
end

