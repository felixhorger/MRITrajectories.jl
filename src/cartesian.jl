
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
	spatial_indices = CartesianIndices(shape)
	indices = Vector{CartesianIndex{N}}(undef, readouts_per_dynamic * num_dynamic)
	for i = 1:num_dynamic
		indices[i:num_dynamic:end] = StatsBase.sample(spatial_indices, readouts_per_dynamic; replace=false)
	end
	return indices
end


"""
	sample_density() = (x,y,z,...) where x, y, z, ... ∈ [-0.5, 0.5]
"""
function variable_density(shape::NTuple{N, Integer}, num::Integer, sample_density::Function; maxiter=100num) where N
	indices = Vector{CartesianIndex{N}}(undef, num)
	i = 1
	j = 0
	upper = CartesianIndex(shape)
	shape_minus_one = shape .- 1
	while i ≤ num
		# Check maximum iterations
		j += 1
		if j == maxiter
			error("Maximum number of iterations exceeded")
			return indices
		end
		# Sample
		sample = sample_density()
		!all(-0.5 .≤ sample .≤ 0.5) && continue
		indices[i] = CartesianIndex(floor.(Int, shape_minus_one .* (sample .+ 0.5)) .+ 1)
		i += 1
	end
	return indices
end

function poisson_disk()
	error("Not implemented")
end

