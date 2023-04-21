
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
	num_samples = num_dynamic * readouts_per_dynamic
	spatial_indices = @views repeat(uniform(shape), ceil(Int, num_samples / prod(shape)))[1:num_samples]
	indices = Vector{CartesianIndex{N}}(undef, readouts_per_dynamic * num_dynamic)
	j = 1
	for i = 1:num_dynamic
		indices[i:num_dynamic:end] = spatial_indices[j:j+readouts_per_dynamic-1]
	end
	return indices
end
#function uniform_dynamic(shape::NTuple{N, Integer}, num_dynamic::Integer, readouts_per_dynamic::Integer) where N
#	spatial_indices = CartesianIndices(shape)
#	indices = Vector{CartesianIndex{N}}(undef, readouts_per_dynamic * num_dynamic)
#	for i = 1:num_dynamic
#		indices[i:num_dynamic:end] = StatsBase.sample(spatial_indices, readouts_per_dynamic; replace=false)
#	end
#	return indices
#end


function tiled(shape::NTuple{N, Integer}, step::NTuple{N, Integer}; random=false) where N
	@assert all(iszero, rem.(shape, step))
	num_tiles = shape .÷ step
	sampling = Vector{CartesianIndex{N}}(undef, prod(shape))
	inner_order = fill(collect(CartesianIndices(step)), num_tiles)
	if random
		shuffle!.(inner_order)
	end
	i = 1
	for I in CartesianIndices(step)
		for J in CartesianIndices(num_tiles)
			sampling[i] = inner_order[J][I] + CartesianIndex{N}((Tuple(J) .- 1) .* step)
			i += 1
		end
	end
	return sampling
end


#function regular_dynamic(shape::NTuple{N, Integer}, num_dynamic::Integer, num_per_k::Integer) where N
#	num_k_per_cycle = num_dynamic ÷ num_per_k
#	step = num_k_per_cycle * num_per_k
#	indices = uniform(shape)
#	num_spatial = prod(shape)
#	num_samples = num_spatial * num_per_k
#	sampling = Vector{CartesianIndex{N}}(undef, num_samples)
#	j = 1
#	i = 1
#	for outer i = 1:step:num_samples
#		(num_samples - i + 1) < step && break
#		k = indices[j:j+num_k_per_cycle-1]
#		j += num_k_per_cycle
#		for l = 0:step-1
#			sampling[i+l] = k[mod(l, num_k_per_cycle)+1]
#		end
#	end
#	k = indices[j:end]
#	length(k) < 5 && @warn "Bad sampling"
#	for l = 0:num_samples-i
#		sampling[i+l] = k[mod(l, num_spatial-j+1)+1]
#	end
#	return sampling
#end

function regular_dynamic(shape::NTuple{N, Integer}, num_dynamic::Integer, num_per_k::Integer) where N
	num_samples = prod(shape) * num_per_k
	sampling = Vector{CartesianIndex{N}}(undef, num_samples)
	to_sample = repeat(uniform(shape), num_per_k)
	j = 1
	for t = 1:num_dynamic
		i = t:num_dynamic:num_samples
		k = j + length(i)
		sampling[i] = to_sample[j:k-1]
		j = k
	end
	return sampling
end




	#num_spatial = prod(shape)
	#indices = vec(CartesianIndices(shape))
	#sampling = Vector{CartesianIndex{N}}(undef, num_spatial * num_per_k)
	#i = 0
	#s = 1
	#sampling .= (CartesianIndex(0,0),)
	#while i < num_spatial
	#	i += 1
	#	@show s, i
	#	calc how many samples can be given, or just min the upper limit that's best I guess
	#	if s+num_dynamic*(num_per_k-1) > num_spatial * num_per_k
	#		break
	#	end
	#	sampling[s:num_dynamic:s+num_dynamic*(num_per_k-1)] .= (indices[i],)
	#	s += 1
	#	if mod1(s, num_dynamic) == 1
	#		s += num_dynamic * (num_per_k-1)
	#	end
	#end

