
"""
"""
# TODO: this belongs into MRITrajectories bc not directly related to reconstruction, but trajectory design
function swap_duplicates_dynamic!(split_indices::AbstractVector{<: AbstractVector{<: CartesianIndex{N}}}; maxiter::Integer=0) where N
	num_dynamic = length(split_indices)
	if maxiter == 0
		maxiter = sum(length.(split_indices))
	end
	for i = 1:num_dynamic
		indices_dynamic = split_indices[i]
		k = zero(CartesianIndex{N})
		j = 1
		while j ≤ length(indices_dynamic)
			m = indices_dynamic[j]
			if k == m
				found = false
				swap_dynamics = randperm!(collect(1:num_dynamic))
				for i_swap in swap_dynamics
					i_swap == i && continue
					indices_dynamic_swap = split_indices[i_swap]
					swap_indices = randperm!(collect(1:length(indices_dynamic_swap)))
					for j_swap in swap_indices
						k_swap = indices_dynamic_swap[j_swap]
						l = searchsortedfirst(indices_dynamic, k_swap)
						l_swap = searchsortedfirst(indices_dynamic_swap, k)
						l ≤ length(indices_dynamic) && indices_dynamic[l] == k_swap && continue # Is k_swap already in indices_dynamic?
						l_swap ≤ length(indices_dynamic_swap) && indices_dynamic_swap[l_swap] == k && continue # Is k already in indices_dynamic_swap?
						# Swap
						insert!(indices_dynamic, l, k_swap)
						insert!(indices_dynamic_swap, l_swap, k)
						deleteat!(indices_dynamic, (l > j) ? j : (j + 1))
						deleteat!(indices_dynamic_swap, (l_swap > j_swap) ? j_swap : (j_swap + 1))
						found = true
						break
					end
					found && break
				end
				if !found
					error("Not implemented, maybe do +- 1 search or reset to a default value, or do a hardcore search")
				end
			else
				j += 1
			end
			k = m
		end
	end
	return split_indices
end


function swap_duplicates_dynamic_nosort!(split_indices::AbstractVector{<: AbstractVector{<: CartesianIndex{N}}}; maxiter::Integer=0) where N
	num_dynamic = length(split_indices)
	if maxiter == 0
		maxiter = sum(length.(split_indices))
	end
	for i = 1:num_dynamic
		indices_dynamic = split_indices[i]
		k = zero(CartesianIndex{N})
		j = 1
		while j ≤ length(indices_dynamic)
			k = indices_dynamic[j]
			if sum(isequal(k), indices_dynamic) > 1
				found = false
				swap_dynamics = randperm!(collect(1:num_dynamic))
				for i_swap in swap_dynamics
					i_swap == i && continue
					indices_dynamic_swap = split_indices[i_swap]
					swap_indices = randperm!(collect(1:length(indices_dynamic_swap)))
					for j_swap in swap_indices
						k_swap = indices_dynamic_swap[j_swap]
						k_swap in indices_dynamic && continue # Is k_swap already in indices_dynamic?
						k in indices_dynamic_swap && continue # Is k already in indices_dynamic_swap?
						# Swap
						indices_dynamic[j], indices_dynamic_swap[j_swap] = k_swap, k
						found = true
						break
					end
					found && break
				end
				if !found
					error("Not implemented, maybe do +- 1 search or reset to a default value, or do a hardcore search")
				end
			else
				j += 1
			end
		end
	end
	return split_indices
end



"""
	floors the floats and adds one
"""
function cartesianise(points::AbstractVector{NTuple{N, Float64}}) where N
	[CartesianIndex(floor.(Int, p) .+ 1) for p in points]
end


function sort_constantL2distance!(split_sampling::AbstractVector{<: AbstractVector{<: CartesianIndex{N}}}, Δ::Real) where N
	@assert length(split_sampling) > 2
	@assert Δ ≥ 0
	num_dynamic = length(split_sampling)
	l = length.(split_sampling)
	for i = 1:sum(l)
		t = mod1(i, num_dynamic)
		j = (i-1) ÷ num_dynamic + 1
		k = Tuple(split_sampling[t][j])
		# Iterate over indices in next dynamic
		r = mod1(t+1, num_dynamic)
		j > l[r] && break
		δ = Inf
		q = 0
		for p = j:l[r]
			δp = abs(Δ - sqrt(sum(abs2, k .- Tuple(split_sampling[r][p]))))
			δp ≥ δ && continue
			δ = δp
			q = p
		end
		# Swap
		split_sampling[r][q], split_sampling[r][j] = split_sampling[r][j], split_sampling[r][q]
	end
	return split_sampling
end

function sort_into_tiles!(split_sampling::AbstractVector{<: AbstractVector{<: CartesianIndex{2}}}, shape::NTuple{2, Integer}, tiles::NTuple{2, Integer})
	# TODO check tile sizes
	tile_size = shape .÷ tiles
	num_dynamic = length(split_sampling)

	spatially_split_sampling = [
		Vector{CartesianIndex{2}}(undef, 0)
		for t = 1:num_dynamic, I in CartesianIndices(tiles)
	]

	for t = 1:num_dynamic, I in split_sampling[t]
		push!(spatially_split_sampling[t, CartesianIndex((Tuple(I) .- 1) .÷ tile_size .+ 1)], I)
	end
	shuffle!.(spatially_split_sampling)
	#@show collect(sum(length.(spatially_split_sampling[t, :, :])) for t = 1:499)
	return collect(sum(length.(spatially_split_sampling[t, I])) for t = 1:499, I in CartesianIndices(tiles))
	

	num_samples = sum(length.(split_sampling))
	sampling = Vector{CartesianIndex{2}}(undef, num_samples)

	direction = ones(Int, 2)
	num_tiles = prod(tiles)
	inside_tile_index = 1
	tile_index = ones(Int, 2)
	i = 1
	while i ≤ num_samples
		t = mod1(i, num_dynamic)
		tile_to_choose_from = spatially_split_sampling[t, tile_index...]
		if length(tile_to_choose_from) > 0
			sampling[i] = pop!(tile_to_choose_from)
			i += 1
		end
		tile_index[1] += direction[1]
		if tile_index[1] > tiles[1] || tile_index[1] == 0
			direction[1] *= -1
			tile_index[1] += direction[1]
			tile_index[2] += direction[2]
			if tile_index[2] > tiles[2] || tile_index[2] == 0
				direction[2] *= -1
				tile_index[2] += direction[2]
			end
		end
	end

	return sampling
end
#function sort_into_tiles!(split_sampling::AbstractVector{<: AbstractVector{<: CartesianIndex{N}}}, shape::NTuple{N, Integer}, tiles::NTuple{N, Integer}) where N
#	# TODO check tile sizes
#	tile_size = shape .÷ tiles
#	num_dynamic = length(split_sampling)
#
#	spatially_split_sampling = [
#		Vector{CartesianIndex{N}}(undef, 0)
#		for t = 1:num_dynamic, I in CartesianIndices(tiles)
#	]
#
#	for t = 1:num_dynamic, I in split_sampling[t]
#		push!(spatially_split_sampling[t, CartesianIndex((Tuple(I) .- 1) .÷ tile_size .+ 1)], I)
#	end
#	shuffle!.(spatially_split_sampling)
#	#@show collect(sum(length.(spatially_split_sampling[t, :, :])) for t = 1:499)
#	#return
#
#	num_samples = sum(length.(split_sampling))
#	sampling = Vector{CartesianIndex{N}}(undef, num_samples)
#
#	direction = 1
#	num_tiles = prod(tiles)
#	tile_index = 1
#	inside_tile_index = 1
#	tile_indices = CartesianIndices(tiles)
#	i = 1
#	while i ≤ num_samples
#		t = mod1(i, num_dynamic)
#		tile_to_choose_from = spatially_split_sampling[t, tile_indices[tile_index]]
#		if length(tile_to_choose_from) > 0
#			sampling[i] = pop!(tile_to_choose_from)
#			i += 1
#		end
#		tile_index += direction
#		if tile_index > num_tiles || tile_index == 0
#			direction *= -1
#			tile_index += direction
#		end
#	end
#
#	return sampling
#end


function linear_order(split_sampling::AbstractVector{<: AbstractVector{<: CartesianIndex{N}}}, shape::NTuple{N, Integer}) where N
	num_dynamic = length(split_sampling)
	split_sampling = copy.(split_sampling)
	num_samples = sum(length.(split_sampling))
	sampling = Vector{CartesianIndex{N}}(undef, num_samples)
	direction = ones(Int, 2)
	index = ones(Int, 2)
	i = 1
	while i ≤ num_samples
		t = mod1(i, num_dynamic)
		#
		mini = Inf
		k = 0
		for (j, I) in enumerate(split_sampling[t])
			d = sum(abs2, Tuple(I) .- index)
			if d < mini
				k = j
			end
		end
		sampling[i] = popat!(split_sampling[t], k)
		i += 1
		#
		index[1] += direction[1]
		if index[1] > shape[1] || index[1] == 0
			direction[1] *= -1
			index[1] += direction[1]
			index[2] += direction[2]
			if index[2] > shape[2] || index[2] == 0
				direction[2] *= -1
				index[2] += direction[2]
			end
		end
	end
	return sampling
end

function sort_increasing(split_sampling::AbstractVector{<: AbstractVector{<: CartesianIndex{N}}}) where N
	num_dynamic = length(split_sampling)
	split_sampling = copy.(split_sampling)
	num_samples = sum(length.(split_sampling))
	sampling = Vector{CartesianIndex{N}}(undef, num_samples)
	i = 1
	while i ≤ num_samples
		t = mod1(i, num_dynamic)
		j = argmin(split_sampling[t])
		sampling[i] = popat!(split_sampling[t], j)
		i += 1
	end
	return sampling
end

