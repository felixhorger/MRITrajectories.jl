
"""
Sample each point in k-space exactly once
"""
function uniform(shape::NTuple{N, Integer}) where N
	indices = vec(collect(CartesianIndices(shape)))
	return shuffle!(indices)
end

"""
	Every point in k-t domain has the same probability
"""
function uniform_dynamic(shape::NTuple{N, Integer}, num_dynamic::Integer, readouts_per_dynamic::Integer) where N
	single_repeat_indices = CartesianIndices(shape)
	indices = Vector{CartesianIndex{N}}(undef, readouts_per_dynamic * num_dynamic)
	for i = 1:num_dynamic
		indices[i:num_dynamic:end] = StatsBase.sample(single_repeat_indices, readouts_per_dynamic; replace=false)
	end
	return indices
end


function poisson_disk()
	error("Not implemented")
end

